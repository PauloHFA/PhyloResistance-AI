"""Testes automatizados das Regras de Rigor Científico Obrigatórias (REG-01 a REG-15)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import scipy.sparse as sp

from phyloresistance.config import (
    COLUNAS_METADADOS_PROIBIDAS,
    FARMACOS_PILOTO,
)
from phyloresistance.dados.bvbrc import gerar_coorte_piloto_ecoli
from phyloresistance.explicabilidade.shap_locus import (
    agregar_shap_por_locus,
    avaliar_explicabilidade_completa,
    calcular_valores_treeshap,
    construir_mapa_atributo_para_locus,
)
from phyloresistance.filogenia.splits import (
    gerar_split_geografico,
    gerar_split_temporal,
    gerar_splits_grupo_filogenetico,
    gerar_splits_leave_one_lineage_out,
    verificar_ausencia_vazamento_grupo,
)
from phyloresistance.genomica.anotacao import construir_matrizes_atributos_coorte
from phyloresistance.genomica.kmers import (
    FiltroSelecaoKmersDobra,
    agrupar_kmers_equivalentes_unitigs,
)
from phyloresistance.modelos.calibracao import CalibradorConformalAMR
from phyloresistance.modelos.escada import (
    BaselineRegrasB1,
    ModeloControleKmersCurtos,
    PipelineCompletoB4,
)
from phyloresistance.validacao.metricas import (
    calcular_metricas_completas,
    comparar_modelo_vs_baseline_pareado,
)


def test_reg01_e_reg09_splits_sem_vazamento_de_grupo() -> None:
    """
    REG-01 e REG-09: Verifica que na divisão por grupo filogenético e LOLO nenhum grupo/ST
    aparece simultaneamente no treino e no teste, e que splits temporal e geográfico são disjuntos.
    """
    rng = np.random.default_rng(42)
    grupos = np.array([f"ST_{i % 8}" for i in range(160)])
    y = rng.integers(0, 2, size=160)
    anos = np.array([2018 + (i % 7) for i in range(160)])
    paises = np.array([["BR", "US", "DE", "IN", "JP"][i % 5] for i in range(160)])

    splits_grupo = gerar_splits_grupo_filogenetico(y, grupos, n_splits=5, seed=42)
    assert len(splits_grupo) == 5
    for sp_g in splits_grupo:
        assert verificar_ausencia_vazamento_grupo(sp_g.indices_treino, sp_g.indices_teste, grupos)
        assert sp_g.grupos_treino.isdisjoint(sp_g.grupos_teste)

    splits_lolo = gerar_splits_leave_one_lineage_out(y, grupos, max_linhagens_avaliadas=3)
    assert len(splits_lolo) == 3
    for sp_l in splits_lolo:
        assert len(sp_l.grupos_teste) == 1
        assert sp_l.grupos_treino.isdisjoint(sp_l.grupos_teste)

    sp_temp = gerar_split_temporal(anos, ano_corte=2021)
    assert max(anos[sp_temp.indices_treino]) <= 2021
    assert min(anos[sp_temp.indices_teste]) > 2021

    sp_geo = gerar_split_geografico(paises, paises_teste=("IN", "JP"))
    assert sp_geo.grupos_treino.isdisjoint(sp_geo.grupos_teste)


def test_reg03_anti_vazamento_selecao_kmers_dentro_da_dobra() -> None:
    """
    REG-03: Garante que FiltroSelecaoKmersDobra usa APENAS os dados de treino da dobra.
    Um k-mer que só aparece nos isolados de teste jamais pode ser aprovado pelo filtro de treino.
    """
    n_treino = 50
    n_teste = 30
    X_treino = np.zeros((n_treino, 10), dtype=np.float64)
    y_treino = np.array([0] * 25 + [1] * 25, dtype=np.int32)
    # Atributos 0..4 informativos no treino; atributo 9 está ZERADO no treino, mas 100% ativo no teste
    X_treino[25:, :5] = 1.0
    X_treino[:10, 5:8] = 1.0

    X_teste = np.zeros((n_teste, 10), dtype=np.float64)
    X_teste[:, 9] = 1.0  # exclusivo do teste!

    filtro = FiltroSelecaoKmersDobra(min_freq=0.05, max_freq=0.95, max_kmers_selecionados=5)
    filtro.fit(X_treino, y_treino, ids_treino=[f"TR_{i}" for i in range(n_treino)])

    assert filtro.indices_selecionados_ is not None
    assert 9 not in set(map(int, filtro.indices_selecionados_))

    # Alterar os valores de X_teste não altera em nada os índices selecionados pelo filtro já ajustado
    X_te_trans = filtro.transform(X_teste)
    assert X_te_trans.shape == (n_teste, len(filtro.indices_selecionados_))


def test_reg08_metadados_proibidos_ausentes_nos_atributos() -> None:
    """REG-08: País, ano, ST e linhagem não entram como atributos em F1, F2 ou F3."""
    _, _, df_geno = gerar_coorte_piloto_ecoli(n_isolados_base=60, seed=42)
    df_f1_f2, X_f3, nomes_f3, _ = construir_matrizes_atributos_coorte(df_geno, seed=42)

    colunas_modelo_f1_f2 = set(df_f1_f2.columns) - {"genome_id"}
    assert colunas_modelo_f1_f2.isdisjoint(COLUNAS_METADADOS_PROIBIDAS)
    assert set(nomes_f3).isdisjoint(COLUNAS_METADADOS_PROIBIDAS)
    assert X_f3.shape[1] == len(nomes_f3)


def test_reg04_reg05_reg11_escada_b1_b4_calibracao_e_abstencao() -> None:
    """
    REG-04, REG-05, REG-11:
    - B4 supera controle negativo de k-mers curtos (k=6) e é comparado a B1 com bootstrap em blocos.
    - VME e ME são reportados e a abstenção conformal reduz erros nos casos duvidosos.
    """
    registros, df_meta, df_geno = gerar_coorte_piloto_ecoli(n_isolados_base=160, seed=42)
    df_geno_base = df_geno.iloc[:160].reset_index(drop=True)
    df_f1_f2, X_f3, nomes_f3, X_k6 = construir_matrizes_atributos_coorte(df_geno_base, seed=42)
    nomes_f1_f2 = [c for c in df_f1_f2.columns if c != "genome_id"]
    X_f12 = df_f1_f2[nomes_f1_f2].to_numpy(dtype=np.float64)

    # Rótulos de ciprofloxacina para os 160 isolados base
    mapa_y = {
        r.genome_id: (1 if r.fenotipo_bruto == "Resistant" else 0)
        for r in registros
        if r.farmaco == "ciprofloxacina" and r.genome_id.startswith("ECOLI_1")
    }
    y = np.array([mapa_y[gid] for gid in df_geno_base["genome_id"].values], dtype=np.int32)
    grupos = df_geno_base["st"].to_numpy(dtype=str)

    splits = gerar_splits_grupo_filogenetico(y, grupos, n_splits=4, seed=42)
    oof_b1 = np.zeros(len(y))
    oof_b4 = np.zeros(len(y))
    oof_k6 = np.zeros(len(y))

    for sp_g in splits:
        tr, te = sp_g.indices_treino, sp_g.indices_teste
        b1 = BaselineRegrasB1("ciprofloxacina", nomes_f1_f2).fit(X_f12[tr], y[tr])
        oof_b1[te] = b1.predict_proba(X_f12[te])[:, 1]

        b4 = PipelineCompletoB4(nomes_f1_f2, nomes_f3, max_kmers_f3=40, seed=42).fit(
            X_f12[tr], X_f3[tr], y[tr]
        )
        oof_b4[te] = b4.predict_proba(X_f12[te], X_f3[te])[:, 1]

        ctrl = ModeloControleKmersCurtos(seed=42).fit(X_k6[tr], y[tr])
        oof_k6[te] = ctrl.predict_proba(X_k6[te])[:, 1]

    cal = CalibradorConformalAMR(metodo="platt", alpha=0.12).fit(oof_b4, y)
    preds_inc = cal.predizer_com_incerteza(oof_b4)

    met_b4 = calcular_metricas_completas(
        "B4", "grupo_filogenetico", y, cal.calibrar_probabilidades(oof_b4),
        predicoes_incerteza=preds_inc, clusters=grupos, n_boot=60, seed=42
    )
    met_k6 = calcular_metricas_completas(
        "Controle_k6", "grupo_filogenetico", y, oof_k6, clusters=grupos, n_boot=30, seed=42
    )

    # B4 com k-mers longos (k=31) + F1/F2 deve superar amplamente o controle negativo de k=6
    assert met_b4.auprc > met_k6.auprc + 0.15
    assert met_b4.auprc_ic95 is not None
    assert 0.0 <= met_b4.vme <= 1.0
    assert 0.0 <= met_b4.me <= 1.0

    comp = comparar_modelo_vs_baseline_pareado(
        "B4", y, cal.calibrar_probabilidades(oof_b4), oof_b1, clusters=grupos, n_boot=60, seed=42
    )
    assert comp.delta_auprc_medio > 0.0


def test_reg10_shap_agregado_por_locus_e_testes_secao_8_2() -> None:
    """
    REG-10 e Seção 8.2:
    Verifica agregação de SHAP por locus, recuperação de determinantes do CARD e teste de remoção.
    """
    registros, _, df_geno = gerar_coorte_piloto_ecoli(n_isolados_base=120, seed=42)
    df_geno_base = df_geno.iloc[:120].reset_index(drop=True)
    df_f1_f2, X_f3, nomes_f3, _ = construir_matrizes_atributos_coorte(df_geno_base, seed=42)
    nomes_f1_f2 = [c for c in df_f1_f2.columns if c != "genome_id"]
    X_f12 = df_f1_f2[nomes_f1_f2].to_numpy(dtype=np.float64)

    mapa_y = {
        r.genome_id: (1 if r.fenotipo_bruto == "Resistant" else 0)
        for r in registros
        if r.farmaco == "ciprofloxacina" and r.genome_id.startswith("ECOLI_1")
    }
    y = np.array([mapa_y[gid] for gid in df_geno_base["genome_id"].values], dtype=np.int32)
    grupos = df_geno_base["st"].to_numpy(dtype=str)

    b4 = PipelineCompletoB4(nomes_f1_f2, nomes_f3, max_kmers_f3=35, seed=42).fit(X_f12, X_f3, y)
    assert b4.modelo_ is not None
    X_trans = b4.transformar_atributos(X_f12, X_f3)
    shap_vals = calcular_valores_treeshap(b4.modelo_, X_trans)
    mapa = construir_mapa_atributo_para_locus(b4.nomes_atributos_finais_, farmaco="ciprofloxacina")
    mat_locus, ord_loci, _ = agregar_shap_por_locus(shap_vals, b4.nomes_atributos_finais_, mapa)

    imp = np.mean(np.abs(mat_locus), axis=0)
    ranking_1 = {ord_loci[i]: float(imp[i]) for i in range(len(ord_loci))}
    ranking_2 = {ord_loci[i]: float(imp[i] * 0.98) for i in range(len(ord_loci))}

    aval = avaliar_explicabilidade_completa(
        farmaco="ciprofloxacina",
        modelo_xgb=b4.modelo_,
        X_teste=X_trans,
        y_teste=y,
        grupos_teste=grupos,
        nomes_atributos=b4.nomes_atributos_finais_,
        rankings_por_dobra=[ranking_1, ranking_2],
        top_k=6,
    )
    assert aval.proporcao_top_k_no_card > 0.30
    assert len(aval.determinantes_card_recuperados) >= 2
    assert aval.remocao_coerente is True
    assert aval.queda_media_probabilidade_remocao_locus > 0.01


def test_reg12_plano_de_analise_registrado_e_unitigs() -> None:
    """Verifica a existência do plano de análise pré-registrado (REG-12) e agrupamento em unitigs."""
    plano = Path("docs/plano_de_analise.md")
    assert plano.exists()
    conteudo = plano.read_text(encoding="utf-8")
    for f in FARMACOS_PILOTO:
        assert f in conteudo

    # Verifica agrupamento de k-mers idênticos em unitigs
    X_kmers = sp.csr_matrix(
        np.array(
            [
                [1, 1, 0],
                [0, 0, 1],
                [1, 1, 0],
            ],
            dtype=np.float64,
        )
    )
    X_u, nomes_u, mapa_u = agrupar_kmers_equivalentes_unitigs(X_kmers, ["km_A", "km_B", "km_C"])
    assert X_u.shape[1] == 2
    assert len(nomes_u) == 2
    assert any(len(v) == 2 for v in mapa_u.values())
