"""Orquestrador ponta a ponta do pipeline PhyloResistance AI (Fases 1 a 7, Multi-Espécie, MIC, B6 e Validação Externa)."""

from __future__ import annotations

import json
import pickle
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import scipy.sparse as sp

from phyloresistance.config import (
    AVISO_USO_PESQUISA,
    ESPECIE_PILOTO,
    ESPECIES_SUPORTADAS,
    FARMACOS_PILOTO,
    FARMACOS_POR_ESPECIE,
    GLOBAL_RANDOM_SEED,
    PONTOS_CORTE_ECOLI,
    VERSOES_BANCOS,
    VERSOES_MODELOS,
    VERSOES_MODELOS_POR_ESPECIE,
)
from phyloresistance.dados.bvbrc import (
    gerar_coorte_especie_extensao,
    gerar_coorte_externa_ncbi_pathogen,
    gerar_coorte_piloto_ecoli,
    gerar_manifesto_dados,
)
from phyloresistance.dados.curadoria import (
    converter_curados_para_dataframe,
    curar_rotulos,
)
from phyloresistance.explicabilidade.shap_locus import (
    AvaliacaoExplicabilidadeQ3,
    agregar_shap_por_locus,
    avaliar_explicabilidade_completa,
    calcular_valores_treeshap,
    construir_mapa_atributo_para_locus,
    extrair_evidencias_isolado,
)
from phyloresistance.figuras import gerar_figuras_svg_relatorio
from phyloresistance.filogenia.distancia import (
    calcular_matriz_distancia_mash,
    desduplicar_genomas_quase_identicos,
    inferir_clusters_filogeneticos,
)
from phyloresistance.filogenia.splits import (
    gerar_split_geografico,
    gerar_split_temporal,
    gerar_splits_aleatorios_vazamento,
    gerar_splits_grupo_filogenetico,
    gerar_splits_leave_one_lineage_out,
)
from phyloresistance.genomica.anotacao import (
    anotar_isolado_fasta,
    construir_matrizes_atributos_coorte,
)
from phyloresistance.genomica.qc import avaliar_qc_coorte
from phyloresistance.modelos.calibracao import CalibradorConformalAMR
from phyloresistance.modelos.escada import (
    BaselineRegrasB1,
    ClassificadorMajoritarioB0,
    ModeloArvoreCuradoB3,
    ModeloControleKmersCurtos,
    ModeloLinearEsparsoKmersB5,
    ModeloLogisticoRegularizadoB2,
    PipelineCompletoB4,
)
from phyloresistance.modelos.extensoes import (
    ModeloMultiSaidaAMR,
    ModeloSequenciaConvolucionalB6,
    RegressorMICLog2,
)
from phyloresistance.validacao.metricas import (
    ComparacaoPareadaBaseline,
    ResultadoMetricasModelo,
    calcular_metricas_completas,
    calcular_metricas_mic_iso20776,
    comparar_modelo_vs_baseline_pareado,
)


@dataclass
class ArtefatoModeloFarmaco:
    """Artefato treinado e calibrado para um par espécie-fármaco, pronto para inferência na API."""

    especie: str
    farmaco: str
    versao_modelo: str
    modelo_b4: PipelineCompletoB4
    calibrador: CalibradorConformalAMR
    nomes_f1_f2: list[str]
    nomes_f3: list[str]
    metricas_validacao_principal: ResultadoMetricasModelo
    comparacao_b1: ComparacaoPareadaBaseline
    avaliacao_xai: AvaliacaoExplicabilidadeQ3
    regressor_mic: RegressorMICLog2 | None = None


class RegistroModelosServico:
    """Gerencia os modelos treinados por espécie e fármaco para o serviço FastAPI."""

    def __init__(self) -> None:
        self.modelos_por_farmaco: dict[str, ArtefatoModeloFarmaco] = {}
        self.modelos_por_especie_farmaco: dict[tuple[str, str], ArtefatoModeloFarmaco] = {}
        self.especies_suportadas: set[str] = set(ESPECIES_SUPORTADAS)
        self.versoes_bancos: dict[str, str] = dict(VERSOES_BANCOS)

    def esta_pronto(self) -> bool:
        return len(self.modelos_por_farmaco) > 0

    def registrar(self, artefato: ArtefatoModeloFarmaco) -> None:
        self.modelos_por_especie_farmaco[(artefato.especie, artefato.farmaco)] = artefato
        self.especies_suportadas.add(artefato.especie)
        if artefato.especie == ESPECIE_PILOTO:
            self.modelos_por_farmaco[artefato.farmaco] = artefato

    def obter_artefato(self, especie: str, farmaco: str) -> ArtefatoModeloFarmaco | None:
        art = self.modelos_por_especie_farmaco.get((especie, farmaco))
        if art is not None:
            return art
        if especie == ESPECIE_PILOTO:
            return self.modelos_por_farmaco.get(farmaco)
        return None

    def listar_todos_artefatos(self) -> list[ArtefatoModeloFarmaco]:
        if self.modelos_por_especie_farmaco:
            return list(self.modelos_por_especie_farmaco.values())
        return list(self.modelos_por_farmaco.values())

    def salvar(self, caminho: str | Path) -> None:
        p = Path(caminho)
        p.parent.mkdir(parents=True, exist_ok=True)
        payload = {
            "por_farmaco": self.modelos_por_farmaco,
            "por_especie_farmaco": self.modelos_por_especie_farmaco,
            "especies_suportadas": self.especies_suportadas,
        }
        with open(p, "wb") as f:
            pickle.dump(payload, f)

    def carregar(self, caminho: str | Path) -> None:
        p = Path(caminho)
        with open(p, "rb") as f:
            dados = pickle.load(f)
        if isinstance(dados, dict) and "por_farmaco" in dados:
            self.modelos_por_farmaco = dados["por_farmaco"]
            self.modelos_por_especie_farmaco = dados.get("por_especie_farmaco", {})
            self.especies_suportadas = set(dados.get("especies_suportadas", ESPECIES_SUPORTADAS))
        elif isinstance(dados, dict):
            self.modelos_por_farmaco = dados
            for f, art in dados.items():
                self.modelos_por_especie_farmaco[(art.especie, f)] = art

    def predizer_isolado_fasta(
        self,
        conteudo_ou_caminho_fasta: str | Path,
        farmacos_solicitados: list[str] | None = None,
        especie: str = ESPECIE_PILOTO,
    ) -> list[dict[str, object]]:
        """
        Realiza a inferência completa para um isolado FASTA:
        1) Anotação genômica F1 + F2 + F3
        2) Predição com modelo B4 + calibração conformal e abstenção + regressão log2(MIC)
        3) Explicação TreeSHAP agregada por locus e cruzada com o CARD/AMRFinderPlus
        """
        atributos_detectados = anotar_isolado_fasta(conteudo_ou_caminho_fasta)
        farmacos_padrao = list(FARMACOS_POR_ESPECIE.get(especie, FARMACOS_PILOTO))
        farmacos = farmacos_solicitados or farmacos_padrao
        predicoes_saida: list[dict[str, object]] = []

        for farmaco in farmacos:
            art = self.obter_artefato(especie, farmaco)
            if art is None:
                predicoes_saida.append(
                    {
                        "farmaco": farmaco,
                        "status": "nao_suportado",
                        "mensagem": f"Par ({especie}, {farmaco}) fora do domínio treinado.",
                    }
                )
                continue

            x_f1_f2 = np.array(
                [[float(atributos_detectados.get(col, 0)) for col in art.nomes_f1_f2]],
                dtype=np.float64,
            )
            x_f3 = sp.csr_matrix(
                np.array(
                    [[float(atributos_detectados.get(col, 0)) for col in art.nomes_f3]],
                    dtype=np.float64,
                )
            )

            prob_bruta = art.modelo_b4.predict_proba(x_f1_f2, x_f3)[:, 1]
            pred_inc = art.calibrador.predizer_com_incerteza(prob_bruta)[0]

            mic_estimado: float | None = None
            if art.regressor_mic is not None:
                mic_arr = art.regressor_mic.predict_mic_mg_l(x_f1_f2, x_f3)
                mic_estimado = float(round(float(mic_arr[0]), 4))

            X_transformado = art.modelo_b4.transformar_atributos(x_f1_f2, x_f3)
            assert art.modelo_b4.modelo_ is not None
            shap_vals = calcular_valores_treeshap(art.modelo_b4.modelo_, X_transformado)
            evidencias = extrair_evidencias_isolado(
                vetor_x_isolado=X_transformado[0],
                vetor_shap_isolado=shap_vals[0],
                nomes_atributos=art.modelo_b4.nomes_atributos_finais_,
                farmaco=farmaco,
                top_k=5,
            )

            predicoes_saida.append(
                {
                    "farmaco": farmaco,
                    "probabilidade_resistencia": pred_inc.probabilidade_resistencia,
                    "classe": pred_inc.classe,
                    "confianca": pred_inc.confianca,
                    "mic_estimado_mg_l": mic_estimado,
                    "evidencias": [
                        {
                            "tipo": ev.tipo,
                            "gene": ev.gene,
                            "alteracao": ev.alteracao,
                            "fonte": ev.fonte,
                            "contribuicao_shap": ev.contribuicao_shap,
                            "status": ev.status,
                        }
                        for ev in evidencias
                    ],
                }
            )
        return predicoes_saida


def _treinar_especie_extensao(
    especie: str,
    registro_servico: RegistroModelosServico,
    n_isolados_base: int = 120,
    seed: int = GLOBAL_RANDOM_SEED,
) -> dict[str, dict[str, object]]:
    """Treina e valida os modelos de extensão da Fase 6 para K. pneumoniae e S. aureus."""
    farmacos = FARMACOS_POR_ESPECIE[especie]
    regs, df_meta, df_geno = gerar_coorte_especie_extensao(
        especie=especie, n_isolados_base=n_isolados_base, seed=seed
    )
    curados, _ = curar_rotulos(regs, especie_alvo=especie, farmacos_alvo=farmacos)
    df_curados = converter_curados_para_dataframe(curados)

    df_f1_f2, X_f3_sparse, nomes_f3, _ = construir_matrizes_atributos_coorte(
        df_geno, seed=seed, especie=especie
    )
    nomes_f1_f2 = [c for c in df_f1_f2.columns if c != "genome_id"]
    X_f12_full = df_f1_f2[nomes_f1_f2].to_numpy(dtype=np.float64)
    mapa_idx = {gid: idx for idx, gid in enumerate(df_geno["genome_id"].values)}

    resumo_especie: dict[str, dict[str, object]] = {}

    for farmaco in farmacos:
        df_sub = df_curados[df_curados["farmaco"] == farmaco].copy()
        idxs = [mapa_idx[gid] for gid in df_sub["genome_id"].values]
        y = df_sub["y"].to_numpy(dtype=np.int32)
        mics = df_sub["mic_mg_l"].to_numpy(dtype=np.float64)
        grupos = df_geno["st"].to_numpy(dtype=str)[idxs]
        X_f12 = X_f12_full[idxs]
        X_f3 = X_f3_sparse[idxs]

        splits = gerar_splits_grupo_filogenetico(y, grupos, n_splits=4, seed=seed)
        oof_b1 = np.zeros(len(y))
        oof_b4 = np.zeros(len(y))

        for sp_g in splits:
            tr, te = sp_g.indices_treino, sp_g.indices_teste
            b1 = BaselineRegrasB1(farmaco=farmaco, nomes_atributos_f1_f2=nomes_f1_f2, especie=especie).fit(
                X_f12[tr], y[tr]
            )
            oof_b1[te] = b1.predict_proba(X_f12[te])[:, 1]

            b4 = PipelineCompletoB4(nomes_f1_f2, nomes_f3, max_kmers_f3=40, seed=seed).fit(
                X_f12[tr], X_f3[tr], y[tr]
            )
            oof_b4[te] = b4.predict_proba(X_f12[te], X_f3[te])[:, 1]

        cal = CalibradorConformalAMR(metodo="platt", alpha=0.12).fit(oof_b4, y)
        probs_cal = cal.calibrar_probabilidades(oof_b4)
        preds_inc = cal.predizer_com_incerteza(oof_b4)

        met_b4 = calcular_metricas_completas(
            nome_modelo="B4_XGBoost_F1_F2_F3",
            estrategia_split="grupo_filogenetico",
            y_true=y,
            probs=probs_cal,
            predicoes_incerteza=preds_inc,
            clusters=grupos,
            n_boot=40,
            seed=seed,
        )
        met_b1 = calcular_metricas_completas(
            nome_modelo="B1_Regras",
            estrategia_split="grupo_filogenetico",
            y_true=y,
            probs=oof_b1,
            n_boot=0,
        )
        comp = comparar_modelo_vs_baseline_pareado(
            "B4_XGBoost_F1_F2_F3", y, probs_cal, oof_b1, clusters=grupos, n_boot=40, seed=seed
        )

        b4_final = PipelineCompletoB4(nomes_f1_f2, nomes_f3, max_kmers_f3=40, seed=seed).fit(
            X_f12, X_f3, y
        )
        reg_mic = RegressorMICLog2(max_kmers_f3=40, seed=seed).fit(X_f12, X_f3, mics)

        assert b4_final.modelo_ is not None
        X_trans = b4_final.transformar_atributos(X_f12, X_f3)
        aval_xai = avaliar_explicabilidade_completa(
            farmaco=farmaco,
            modelo_xgb=b4_final.modelo_,
            X_teste=X_trans,
            y_teste=y,
            grupos_teste=grupos,
            nomes_atributos=b4_final.nomes_atributos_finais_,
            rankings_por_dobra=[],
            top_k=4,
        )

        versao_m = VERSOES_MODELOS_POR_ESPECIE[especie][farmaco]
        art = ArtefatoModeloFarmaco(
            especie=especie,
            farmaco=farmaco,
            versao_modelo=versao_m,
            modelo_b4=b4_final,
            calibrador=cal,
            nomes_f1_f2=nomes_f1_f2,
            nomes_f3=nomes_f3,
            metricas_validacao_principal=met_b4,
            comparacao_b1=comp,
            avaliacao_xai=aval_xai,
            regressor_mic=reg_mic,
        )
        registro_servico.registrar(art)

        resumo_especie[farmaco] = {
            "versao_modelo": versao_m,
            "total_isolados": len(y),
            "auprc_b1": met_b1.auprc,
            "auprc_b4": met_b4.auprc,
            "vme_b4": met_b4.vme,
            "me_b4": met_b4.me,
            "delta_auprc_vs_b1": comp.delta_auprc_medio,
            "supera_b1": comp.supera_b1_com_ic95,
            "principais_loci_shap": [loc["locus"] for loc in aval_xai.principais_loci_globais[:3]],
        }

    return resumo_especie


def executar_pipeline_completo(
    output_dir: str | Path = "data/processed",
    reports_dir: str | Path = "reports",
    n_isolados_base: int = 360,
    n_boot: int = 200,
    seed: int = GLOBAL_RANDOM_SEED,
) -> tuple[RegistroModelosServico, dict[str, object]]:
    """
    Executa de ponta a ponta as Etapas 1 a 8 da arquitetura do PhyloResistance AI (Fases 1 a 7):
    [1] Ingestão e curadoria padronizada de rótulos (mais análise de sensibilidade de 'I')
    [2] QC dos genomas + confirmação de espécie + desduplicação clonal (REG-02, REG-07)
    [3] Estrutura filogenética (distâncias Mash -> clusters / ST)
    [4] Engenharia de atributos: F1 (genes AMR), F2 (mutações pontuais), F3 (k-mers k=31 esparsos) e Controle k=6
    [5] Divisão por grupo filogenético (principal), LOLO, Temporal, Geográfica e Aleatória
    [6] Escada de modelos (B0, B1, B2, B3, B4, B5, B6 CNN1D e Controle k=6)
    [7] Calibração probabilística + incerteza / abstenção conformal (Q4, REG-11)
    [8] Explicabilidade TreeSHAP agregada por locus + 4 testes de avaliação da Seção 8.2 (Q3, REG-10)
    [9] Fase 6: Regressão log2(MIC) ISO 20776-2, Multi-Saída, Multi-Espécie (K. pneumoniae, S. aureus) e Validação Externa (NCBI Pathogen Detection)
    [10] Fase 7: Figuras SVG vetoriais e Relatórios Quarto / Markdown
    """
    dir_proc = Path(output_dir)
    dir_rep = Path(reports_dir)
    dir_proc.mkdir(parents=True, exist_ok=True)
    dir_rep.mkdir(parents=True, exist_ok=True)

    # --- ETAPA 1: Ingestão e Curadoria de Rótulos ---
    registros_brutos, df_genomas_meta, df_genotipos_verdade = gerar_coorte_piloto_ecoli(
        n_isolados_base=n_isolados_base, seed=seed
    )
    curados_principal, funil = curar_rotulos(
        registros_brutos,
        especie_alvo=ESPECIE_PILOTO,
        farmacos_alvo=FARMACOS_PILOTO,
        modo_intermediario="excluir",
    )
    curados_sens_i_como_s, _ = curar_rotulos(
        registros_brutos,
        especie_alvo=ESPECIE_PILOTO,
        farmacos_alvo=FARMACOS_PILOTO,
        modo_intermediario="agrupar_s",
    )
    df_curados = converter_curados_para_dataframe(curados_principal)

    # --- ETAPA 2: QC dos Genomas + Confirmação de Espécie ---
    df_qc, _ = avaliar_qc_coorte(df_genomas_meta)
    ids_reprovados_qc = set(df_qc.loc[~df_qc["aprovado"], "genome_id"].values)
    funil.removidos_qc_genoma = int(df_curados["genome_id"].isin(ids_reprovados_qc).sum())

    df_aprovados_qc = df_qc[df_qc["aprovado"]].reset_index(drop=True)
    ids_aprovados_qc = list(df_aprovados_qc["genome_id"].values)

    df_genotipos_qc = (
        df_genotipos_verdade.set_index("genome_id")
        .loc[ids_aprovados_qc]
        .reset_index()
    )

    # --- ETAPA 3 e 4 (Pré-Split): Atributos Base, Distâncias Mash e Desduplicação (REG-02) ---
    df_f1_f2_bruto, X_f3_bruto, nomes_f3, X_k6_bruto = construir_matrizes_atributos_coorte(
        df_genotipos_qc, seed=seed
    )
    X_comb_dist = sp.hstack(
        [sp.csr_matrix(df_f1_f2_bruto.drop(columns=["genome_id"]).to_numpy(dtype=np.float64)), X_f3_bruto]
    ).tocsr()
    matriz_dist_mash = calcular_matriz_distancia_mash(X_comb_dist, k=31)

    indices_unicos, relatorio_dedup = desduplicar_genomas_quase_identicos(
        genome_ids=ids_aprovados_qc,
        matriz_distancia=matriz_dist_mash,
        limiar_distancia=0.0005,
    )
    ids_removidos_dedup = set(relatorio_dedup.ids_removidos)
    funil.removidos_desduplicacao_clonal = int(df_curados["genome_id"].isin(ids_removidos_dedup).sum())

    df_meta_unico = df_aprovados_qc.iloc[indices_unicos].reset_index(drop=True)
    df_f1_f2 = df_f1_f2_bruto.iloc[indices_unicos].reset_index(drop=True)
    X_f3_sparse = X_f3_bruto[indices_unicos]
    X_k6_controle = X_k6_bruto[indices_unicos]
    matriz_dist_unica = matriz_dist_mash[np.ix_(indices_unicos, indices_unicos)]

    clusters_filogeneticos = inferir_clusters_filogeneticos(
        genome_ids=list(df_meta_unico["genome_id"].values),
        matriz_distancia=matriz_dist_unica,
        sts=list(df_meta_unico["st"].values),
    )
    df_meta_unico["cluster_filogenetico"] = clusters_filogeneticos

    ids_finais_validos = set(df_meta_unico["genome_id"].values)
    df_curados_final = df_curados[df_curados["genome_id"].isin(ids_finais_validos)].reset_index(drop=True)
    funil.total_aprovados_pos_qc_dedup = len(df_curados_final)

    caminho_rotulos = dir_proc / "curadoria_rotulos.parquet"
    caminho_qc = dir_proc / "qc_genomas.parquet"
    caminho_clusters = dir_proc / "clusters_filogeneticos.parquet"
    caminho_f1_f2 = dir_proc / "atributos_f1_f2.parquet"
    caminho_f3 = dir_proc / "atributos_f3_kmers.npz"
    caminho_funil = dir_proc / "funil_exclusao.json"
    caminho_dedup = dir_proc / "relatorio_desduplicacao.json"

    df_curados_final.to_parquet(caminho_rotulos, index=False)
    df_qc.to_parquet(caminho_qc, index=False)
    df_meta_unico[["genome_id", "st", "cluster_filogenetico", "pais", "ano"]].to_parquet(
        caminho_clusters, index=False
    )
    df_f1_f2.to_parquet(caminho_f1_f2, index=False)
    sp.save_npz(caminho_f3, X_f3_sparse)
    caminho_funil.write_text(json.dumps(funil.model_dump(), indent=2, ensure_ascii=False), encoding="utf-8")
    caminho_dedup.write_text(
        json.dumps(relatorio_dedup.model_dump(), indent=2, ensure_ascii=False), encoding="utf-8"
    )

    # --- ETAPAS 5 a 8: Modelagem, Validação, Calibração, XAI e Regressão de MIC por Fármaco ---
    nomes_f1_f2 = [c for c in df_f1_f2.columns if c != "genome_id"]
    X_f1_f2_full = df_f1_f2[nomes_f1_f2].to_numpy(dtype=np.float64)
    mapa_id_para_idx = {gid: idx for idx, gid in enumerate(df_meta_unico["genome_id"].values)}

    registro_servico = RegistroModelosServico()
    resultados_por_farmaco: dict[str, object] = {}
    resultados_mic_iso: dict[str, dict[str, float | str]] = {}

    for farmaco in FARMACOS_PILOTO:
        df_sub = df_curados_final[df_curados_final["farmaco"] == farmaco].copy()
        indices_alinhados = [mapa_id_para_idx[gid] for gid in df_sub["genome_id"].values]
        y = df_sub["y"].to_numpy(dtype=np.int32)
        mics_reais = df_sub["mic_mg_l"].to_numpy(dtype=np.float64)
        ids_isolados = list(df_sub["genome_id"].values)
        grupos = np.array([clusters_filogeneticos[i] for i in indices_alinhados], dtype=str)
        anos = df_meta_unico["ano"].to_numpy(dtype=np.int32)[indices_alinhados]
        paises = df_meta_unico["pais"].to_numpy(dtype=str)[indices_alinhados]

        X_f12 = X_f1_f2_full[indices_alinhados]
        X_f3 = X_f3_sparse[indices_alinhados]
        X_k6 = X_k6_controle[indices_alinhados]

        splits_grupo = gerar_splits_grupo_filogenetico(y, grupos, n_splits=5, seed=seed)
        n = len(y)
        oof_probs: dict[str, np.ndarray] = {
            "B0_Classe_Majoritaria": np.zeros(n, dtype=np.float64),
            "B1_Regras_CARD_AMRFinder": np.zeros(n, dtype=np.float64),
            "B2_LogReg_L1_F1_F2": np.zeros(n, dtype=np.float64),
            "B3_XGBoost_F1_F2": np.zeros(n, dtype=np.float64),
            "B4_XGBoost_F1_F2_F3": np.zeros(n, dtype=np.float64),
            "B5_Linear_Esparso_F3": np.zeros(n, dtype=np.float64),
            "B6_CNN1D_Sequencia": np.zeros(n, dtype=np.float64),
            "Controle_Kmers_Curtos_k6": np.zeros(n, dtype=np.float64),
        }
        oof_probs_b4_brutas = np.zeros(n, dtype=np.float64)
        oof_log2_mic = np.zeros(n, dtype=np.float64)
        rankings_shap_por_dobra: list[dict[str, float]] = []

        for split in splits_grupo:
            tr, te = split.indices_treino, split.indices_teste
            ids_tr = [ids_isolados[int(i)] for i in tr]

            m_b0 = ClassificadorMajoritarioB0().fit(X_f12[tr], y[tr])
            oof_probs["B0_Classe_Majoritaria"][te] = m_b0.predict_proba(X_f12[te])[:, 1]

            m_b1 = BaselineRegrasB1(farmaco=farmaco, nomes_atributos_f1_f2=nomes_f1_f2).fit(X_f12[tr], y[tr])
            oof_probs["B1_Regras_CARD_AMRFinder"][te] = m_b1.predict_proba(X_f12[te])[:, 1]

            m_b2 = ModeloLogisticoRegularizadoB2(seed=seed).fit(X_f12[tr], y[tr])
            oof_probs["B2_LogReg_L1_F1_F2"][te] = m_b2.predict_proba(X_f12[te])[:, 1]

            m_b3 = ModeloArvoreCuradoB3(seed=seed).fit(X_f12[tr], y[tr])
            oof_probs["B3_XGBoost_F1_F2"][te] = m_b3.predict_proba(X_f12[te])[:, 1]

            m_b4 = PipelineCompletoB4(
                nomes_f1_f2=nomes_f1_f2,
                nomes_f3=nomes_f3,
                max_kmers_f3=50,
                seed=seed,
            ).fit(X_f12[tr], X_f3[tr], y[tr], ids_treino=ids_tr)
            probs_te_b4 = m_b4.predict_proba(X_f12[te], X_f3[te])[:, 1]
            oof_probs_b4_brutas[te] = probs_te_b4
            oof_probs["B4_XGBoost_F1_F2_F3"][te] = probs_te_b4

            assert m_b4.modelo_ is not None
            X_te_b4_trans = m_b4.transformar_atributos(X_f12[te], X_f3[te])
            shap_dobra = calcular_valores_treeshap(m_b4.modelo_, X_te_b4_trans)
            mapa_dobra = construir_mapa_atributo_para_locus(
                m_b4.nomes_atributos_finais_, farmaco=farmaco
            )
            mat_loc_dobra, ord_loc_dobra, _ = agregar_shap_por_locus(
                shap_dobra, m_b4.nomes_atributos_finais_, mapa_dobra
            )
            imp_dobra = np.mean(np.abs(mat_loc_dobra), axis=0)
            rankings_shap_por_dobra.append(
                {ord_loc_dobra[i]: float(imp_dobra[i]) for i in range(len(ord_loc_dobra))}
            )

            m_b5 = ModeloLinearEsparsoKmersB5(max_kmers_f3=60, seed=seed).fit(
                X_f3[tr], y[tr], ids_treino=ids_tr
            )
            oof_probs["B5_Linear_Esparso_F3"][te] = m_b5.predict_proba(X_f3[te])[:, 1]

            m_b6 = ModeloSequenciaConvolucionalB6(seed=seed).fit(X_f12[tr], X_f3[tr], y[tr])
            oof_probs["B6_CNN1D_Sequencia"][te] = m_b6.predict_proba(X_f12[te], X_f3[te])[:, 1]

            m_ctrl = ModeloControleKmersCurtos(seed=seed).fit(X_k6[tr], y[tr])
            oof_probs["Controle_Kmers_Curtos_k6"][te] = m_ctrl.predict_proba(X_k6[te])[:, 1]

            reg_m = RegressorMICLog2(max_kmers_f3=50, seed=seed).fit(X_f12[tr], X_f3[tr], mics_reais[tr])
            oof_log2_mic[te] = reg_m.predict_log2_mic(X_f12[te], X_f3[te])

        # Avaliação de Regressão log2(MIC) ISO 20776-2 (Fase 6)
        aval_mic = calcular_metricas_mic_iso20776(
            farmaco=farmaco,
            mic_real_mg_l=mics_reais,
            log2_mic_pred=oof_log2_mic,
            ponto_corte=PONTOS_CORTE_ECOLI.get(farmaco),
        )
        resultados_mic_iso[farmaco] = aval_mic.model_dump()

        # Calibração Out-of-Fold e Predição Conformal (Q4, REG-11)
        calibrador_oof = CalibradorConformalAMR(metodo="platt", alpha=0.12)
        calibrador_oof.fit(oof_probs_b4_brutas, y)
        probs_b4_calibradas = calibrador_oof.calibrar_probabilidades(oof_probs_b4_brutas)
        oof_probs["B4_XGBoost_F1_F2_F3"] = probs_b4_calibradas
        preds_incerteza_b4 = calibrador_oof.predizer_com_incerteza(oof_probs_b4_brutas)

        metricas_escada: dict[str, dict[str, object]] = {}
        for nome_mod, p_vec in oof_probs.items():
            res_m = calcular_metricas_completas(
                nome_modelo=nome_mod,
                estrategia_split="grupo_filogenetico",
                y_true=y,
                probs=p_vec,
                predicoes_incerteza=preds_incerteza_b4 if nome_mod == "B4_XGBoost_F1_F2_F3" else None,
                clusters=grupos,
                n_boot=n_boot if nome_mod in {"B1_Regras_CARD_AMRFinder", "B4_XGBoost_F1_F2_F3"} else 40,
                seed=seed,
            )
            metricas_escada[nome_mod] = res_m.model_dump()

        comp_b3_vs_b1 = comparar_modelo_vs_baseline_pareado(
            nome_modelo="B3_XGBoost_F1_F2",
            y_true=y,
            probs_modelo=oof_probs["B3_XGBoost_F1_F2"],
            probs_b1=oof_probs["B1_Regras_CARD_AMRFinder"],
            clusters=grupos,
            n_boot=n_boot,
            seed=seed,
        )
        comp_b4_vs_b1 = comparar_modelo_vs_baseline_pareado(
            nome_modelo="B4_XGBoost_F1_F2_F3",
            y_true=y,
            probs_modelo=oof_probs["B4_XGBoost_F1_F2_F3"],
            probs_b1=oof_probs["B1_Regras_CARD_AMRFinder"],
            clusters=grupos,
            n_boot=n_boot,
            seed=seed,
        )

        # Validações de Generalização (Q2 e REG-09)
        split_temp = gerar_split_temporal(anos, ano_corte=2022)
        b4_temp = PipelineCompletoB4(nomes_f1_f2=nomes_f1_f2, nomes_f3=nomes_f3, seed=seed).fit(
            X_f12[split_temp.indices_treino],
            X_f3[split_temp.indices_treino],
            y[split_temp.indices_treino],
        )
        prob_temp = b4_temp.predict_proba(
            X_f12[split_temp.indices_teste], X_f3[split_temp.indices_teste]
        )[:, 1]
        met_temporal = calcular_metricas_completas(
            nome_modelo="B4_XGBoost_F1_F2_F3",
            estrategia_split="temporal",
            y_true=y[split_temp.indices_teste],
            probs=prob_temp,
            n_boot=0,
        )

        split_geo = gerar_split_geografico(paises, paises_teste=("IN", "JP", "FR"))
        b4_geo = PipelineCompletoB4(nomes_f1_f2=nomes_f1_f2, nomes_f3=nomes_f3, seed=seed).fit(
            X_f12[split_geo.indices_treino],
            X_f3[split_geo.indices_treino],
            y[split_geo.indices_treino],
        )
        prob_geo = b4_geo.predict_proba(
            X_f12[split_geo.indices_teste], X_f3[split_geo.indices_teste]
        )[:, 1]
        met_geografica = calcular_metricas_completas(
            nome_modelo="B4_XGBoost_F1_F2_F3",
            estrategia_split="geografica",
            y_true=y[split_geo.indices_teste],
            probs=prob_geo,
            n_boot=0,
        )

        splits_lolo = gerar_splits_leave_one_lineage_out(y, grupos, max_linhagens_avaliadas=4)
        y_lolo_all: list[np.ndarray] = []
        p_lolo_all: list[np.ndarray] = []
        for sp_l in splits_lolo:
            b4_lolo = PipelineCompletoB4(nomes_f1_f2=nomes_f1_f2, nomes_f3=nomes_f3, seed=seed).fit(
                X_f12[sp_l.indices_treino],
                X_f3[sp_l.indices_treino],
                y[sp_l.indices_treino],
            )
            p_lolo = b4_lolo.predict_proba(
                X_f12[sp_l.indices_teste], X_f3[sp_l.indices_teste]
            )[:, 1]
            y_lolo_all.append(y[sp_l.indices_teste])
            p_lolo_all.append(p_lolo)
        met_lolo = calcular_metricas_completas(
            nome_modelo="B4_XGBoost_F1_F2_F3",
            estrategia_split="leave_one_lineage_out",
            y_true=np.concatenate(y_lolo_all),
            probs=np.concatenate(p_lolo_all),
            n_boot=0,
        )

        splits_aleat = gerar_splits_aleatorios_vazamento(y, grupos, n_splits=5, seed=seed)
        oof_b5_aleat = np.zeros(n, dtype=np.float64)
        for sp_a in splits_aleat:
            b5_a = ModeloLinearEsparsoKmersB5(max_kmers_f3=60, seed=seed).fit(
                X_f3[sp_a.indices_treino], y[sp_a.indices_treino]
            )
            oof_b5_aleat[sp_a.indices_teste] = b5_a.predict_proba(X_f3[sp_a.indices_teste])[:, 1]
        met_vazamento_aleatorio = calcular_metricas_completas(
            nome_modelo="B5_Linear_Esparso_F3 (Split Aleatório com Vazamento)",
            estrategia_split="aleatoria_vazamento",
            y_true=y,
            probs=oof_b5_aleat,
            n_boot=0,
        )

        modelo_b4_final = PipelineCompletoB4(
            nomes_f1_f2=nomes_f1_f2,
            nomes_f3=nomes_f3,
            max_kmers_f3=50,
            seed=seed,
        ).fit(X_f12, X_f3, y, ids_treino=ids_isolados)
        reg_mic_final = RegressorMICLog2(max_kmers_f3=50, seed=seed).fit(X_f12, X_f3, mics_reais)

        assert modelo_b4_final.modelo_ is not None
        X_full_trans = modelo_b4_final.transformar_atributos(X_f12, X_f3)
        avaliacao_xai = avaliar_explicabilidade_completa(
            farmaco=farmaco,
            modelo_xgb=modelo_b4_final.modelo_,
            X_teste=X_full_trans,
            y_teste=y,
            grupos_teste=grupos,
            nomes_atributos=modelo_b4_final.nomes_atributos_finais_,
            rankings_por_dobra=rankings_shap_por_dobra,
            top_k=6,
        )

        metricas_b4_obj = ResultadoMetricasModelo(**metricas_escada["B4_XGBoost_F1_F2_F3"])
        artefato = ArtefatoModeloFarmaco(
            especie=ESPECIE_PILOTO,
            farmaco=farmaco,
            versao_modelo=VERSOES_MODELOS[farmaco],
            modelo_b4=modelo_b4_final,
            calibrador=calibrador_oof,
            nomes_f1_f2=nomes_f1_f2,
            nomes_f3=nomes_f3,
            metricas_validacao_principal=metricas_b4_obj,
            comparacao_b1=comp_b4_vs_b1,
            avaliacao_xai=avaliacao_xai,
            regressor_mic=reg_mic_final,
        )
        registro_servico.registrar(artefato)

        resultados_por_farmaco[farmaco] = {
            "versao_modelo": VERSOES_MODELOS[farmaco],
            "total_isolados": int(n),
            "prevalencia_resistencia": round(float(np.mean(y)), 4),
            "escada_modelos_grupo_filogenetico": metricas_escada,
            "comparacao_pareada_vs_b1": {
                "B3_vs_B1": comp_b3_vs_b1.model_dump(),
                "B4_vs_B1": comp_b4_vs_b1.model_dump(),
            },
            "generalizacao_q2": {
                "grupo_filogenetico": metricas_escada["B4_XGBoost_F1_F2_F3"],
                "leave_one_lineage_out": met_lolo.model_dump(),
                "temporal": met_temporal.model_dump(),
                "geografica": met_geografica.model_dump(),
                "demonstracao_vazamento_aleatorio_b5": met_vazamento_aleatorio.model_dump(),
            },
            "explicabilidade_q3": avaliacao_xai.model_dump(),
            "regressao_mic_iso20776": aval_mic.model_dump(),
        }

    # --- FASE 6: Extensões (Multi-Espécie, Multi-Saída e Validação Externa NCBI Pathogen Detection) ---
    # 1) Multi-Espécie: Klebsiella pneumoniae (Fase 2) e Staphylococcus aureus (Fase 3)
    res_kpneu = _treinar_especie_extensao(
        "Klebsiella pneumoniae", registro_servico, n_isolados_base=max(90, n_isolados_base // 3), seed=seed
    )
    res_saureus = _treinar_especie_extensao(
        "Staphylococcus aureus", registro_servico, n_isolados_base=max(90, n_isolados_base // 3), seed=seed
    )

    # 2) Validação Externa Independente (NCBI Pathogen Detection AST — Critério de Aceite da Fase 6)
    regs_ext, _, df_geno_ext = gerar_coorte_externa_ncbi_pathogen(
        especie=ESPECIE_PILOTO, n_isolados=100, seed=seed + 1000
    )
    curados_ext, _ = curar_rotulos(regs_ext, especie_alvo=ESPECIE_PILOTO, farmacos_alvo=FARMACOS_PILOTO)
    df_cur_ext = converter_curados_para_dataframe(curados_ext)
    df_f12_ext, X_f3_ext, _, _ = construir_matrizes_atributos_coorte(
        df_geno_ext, seed=seed + 1000, especie=ESPECIE_PILOTO
    )
    X_f12_ext_arr = df_f12_ext[nomes_f1_f2].to_numpy(dtype=np.float64)
    mapa_ext = {gid: idx for idx, gid in enumerate(df_geno_ext["genome_id"].values)}

    validacao_externa_ncbi: dict[str, dict[str, float]] = {}
    for farmaco in FARMACOS_PILOTO:
        df_sub_e = df_cur_ext[df_cur_ext["farmaco"] == farmaco]
        idxs_e = [mapa_ext[gid] for gid in df_sub_e["genome_id"].values]
        y_ext = df_sub_e["y"].to_numpy(dtype=np.int32)
        art_f = registro_servico.modelos_por_farmaco[farmaco]

        probs_b4_ext = art_f.calibrador.calibrar_probabilidades(
            art_f.modelo_b4.predict_proba(X_f12_ext_arr[idxs_e], X_f3_ext[idxs_e])[:, 1]
        )
        b1_ext = BaselineRegrasB1(farmaco=farmaco, nomes_atributos_f1_f2=nomes_f1_f2).fit(
            X_f12_ext_arr[idxs_e], y_ext
        )
        probs_b1_ext = b1_ext.predict_proba(X_f12_ext_arr[idxs_e])[:, 1]

        met_b4_ext = calcular_metricas_completas(
            "B4", "externa_ncbi", y_ext, probs_b4_ext, n_boot=0
        )
        met_b1_ext = calcular_metricas_completas(
            "B1", "externa_ncbi", y_ext, probs_b1_ext, n_boot=0
        )
        validacao_externa_ncbi[farmaco] = {
            "auprc_b1_externo": met_b1_ext.auprc,
            "auprc_b4_externo": met_b4_ext.auprc,
            "ganho_auprc_b4_sobre_b1": round(met_b4_ext.auprc - met_b1_ext.auprc, 4),
            "vme_b4_externo": met_b4_ext.vme,
            "me_b4_externo": met_b4_ext.me,
        }

    # 3) Modelo Multi-Saída Opcional (Seção 2)
    tabela_pivo = df_curados_final.pivot(index="genome_id", columns="farmaco", values="y").dropna()
    ids_pivo = list(tabela_pivo.index)
    idxs_pivo = [mapa_id_para_idx[gid] for gid in ids_pivo]
    Y_multi = tabela_pivo[list(FARMACOS_PILOTO)].to_numpy(dtype=np.int32)
    X_multi = X_f1_f2_full[idxs_pivo]
    mod_multi = ModeloMultiSaidaAMR(seed=seed).fit(X_multi, Y_multi)
    probs_multi = mod_multi.predict_proba_matriz(X_multi)
    auprc_medio_multisaida = float(
        np.mean(
            [
                calcular_metricas_completas("MultiSaida", "treino", Y_multi[:, j], probs_multi[:, j], n_boot=0).auprc
                for j in range(len(FARMACOS_PILOTO))
            ]
        )
    )

    extensoes_fase6: dict[str, object] = {
        "multiespecie": {
            "Klebsiella pneumoniae": res_kpneu,
            "Staphylococcus aureus": res_saureus,
        },
        "regressao_mic_iso20776": resultados_mic_iso,
        "validacao_externa_ncbi_pathogen": validacao_externa_ncbi,
        "modelo_multisaida_auprc_medio": round(auprc_medio_multisaida, 4),
    }

    caminho_modelos = dir_proc / "modelos_servico.pkl"
    registro_servico.salvar(caminho_modelos)

    manifesto = gerar_manifesto_dados(
        [
            caminho_rotulos,
            caminho_qc,
            caminho_clusters,
            caminho_f1_f2,
            caminho_f3,
            caminho_funil,
            caminho_dedup,
            caminho_modelos,
        ]
    )
    caminho_manifesto = dir_proc / "manifesto_versoes_hashes.json"
    caminho_manifesto.write_text(
        json.dumps(manifesto.model_dump(), indent=2, ensure_ascii=False), encoding="utf-8"
    )

    relatorio_consolidado: dict[str, object] = {
        "projeto": "PhyloResistance AI",
        "especie_piloto": ESPECIE_PILOTO,
        "especies_suportadas": list(ESPECIES_SUPORTADAS),
        "farmacos_piloto": list(FARMACOS_PILOTO),
        "versoes_bancos": dict(VERSOES_BANCOS),
        "funil_exclusao": funil.model_dump(),
        "analise_sensibilidade_intermediario": {
            "total_excluindo_I": len(curados_principal),
            "total_agrupando_I_com_S": len(curados_sens_i_como_s),
        },
        "desduplicacao": relatorio_dedup.model_dump(),
        "resultados_por_farmaco": resultados_por_farmaco,
        "extensoes_fase6": extensoes_fase6,
        "manifesto_dados": manifesto.model_dump(),
        "aviso": AVISO_USO_PESQUISA,
    }

    # Gerar Figuras SVG Vetoriais da Fase 7
    caminhos_svg = gerar_figuras_svg_relatorio(
        relatorio_consolidado, diretorio_figuras=dir_rep / "figuras"
    )
    relatorio_consolidado["figuras_svg_geradas"] = caminhos_svg

    caminho_rel_json = dir_rep / "relatorio_validacao.json"
    caminho_rel_json.write_text(
        json.dumps(relatorio_consolidado, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    caminho_rel_md = dir_rep / "relatorio_cientifico.md"
    caminho_rel_md.write_text(
        gerar_relatorio_cientifico_markdown(relatorio_consolidado), encoding="utf-8"
    )

    return registro_servico, relatorio_consolidado


def gerar_relatorio_cientifico_markdown(rel: dict[str, object]) -> str:
    """Gera o relatório científico final (formato Quarto / Markdown) respondendo Q1 a Q4 + Extensões Fase 6."""
    funil = rel["funil_exclusao"]  # type: ignore[index]
    dedup = rel["desduplicacao"]  # type: ignore[index]
    res_farmacos = rel["resultados_por_farmaco"]  # type: ignore[index]
    ext6 = rel.get("extensoes_fase6", {})

    secoes_farmacos: list[str] = []
    for farmaco, dados in res_farmacos.items():  # type: ignore[union-attr]
        escada = dados["escada_modelos_grupo_filogenetico"]
        comp_b4 = dados["comparacao_pareada_vs_b1"]["B4_vs_B1"]
        gen_q2 = dados["generalizacao_q2"]
        xai_q3 = dados["explicabilidade_q3"]
        mic_iso = dados.get("regressao_mic_iso20776", {})

        linhas_escada = []
        for nome_m, m in escada.items():
            ic_str = (
                f"[{m['auprc_ic95']['ic95_inferior']:.3f}, {m['auprc_ic95']['ic95_superior']:.3f}]"
                if m.get("auprc_ic95")
                else "-"
            )
            linhas_escada.append(
                f"| `{nome_m}` | {m['auprc']:.4f} {ic_str} | {m['auroc']:.4f} | "
                f"{m['sensibilidade']:.4f} | {m['especificidade']:.4f} | "
                f"**{m['vme']:.4f}** | **{m['me']:.4f}** | {m['brier_score']:.4f} |"
            )

        linhas_loci = []
        for loc in xai_q3["principais_loci_globais"]:
            linhas_loci.append(
                f"| `{loc['locus']}` | {loc['tipo']} | **{loc['status']}** | {loc['fonte']} | "
                f"{loc['importancia_shap_media']:.4f} | {loc['mecanismo']} |"
            )

        m_b4 = escada["B4_XGBoost_F1_F2_F3"]
        secoes_farmacos.append(
            f"""### Fármaco: `{farmaco}` (`{dados['versao_modelo']}`)

- **Isolados únicos pós-QC/desduplicação**: `{dados['total_isolados']}` (Prevalência de Resistência: `{dados['prevalencia_resistencia']*100:.1f}%`)
- **Conclusão Q1 (ML vs Baseline B1 em grupos filogenéticos não vistos)**: {comp_b4['conclusao_q1']}
- **Regressão Quantitativa de $\\log_2(\\text{{MIC}})$ (ISO 20776-2)**: Concordância Essencial ($\\pm 1$ diluição) = `{mic_iso.get('concordancia_essencial_1_diluicao', 0.0)*100:.1f}%` | Concordância Categórica = `{mic_iso.get('concordancia_categorica', 0.0)*100:.1f}%` | RMSE $\\log_2(\\text{{MIC}})$ = `{mic_iso.get('rmse_log2_mic', 0.0):.3f}`

#### Escada de Modelos (Validação Principal: Grupo Filogenético — REG-01, REG-04, REG-05)

| Modelo | AUPRC (IC 95% Bloco) | AUROC | Sensibilidade | Especificidade | VME (Erro Muito Grave) | ME (Erro Grave) | Brier Score |
|---|---|---|---|---|---|---|---|
{chr(10).join(linhas_escada)}

#### Generalização entre Estratégias de Split (Q2 e REG-09)

| Estratégia de Validação | AUPRC | AUROC | VME | ME | Observação |
|---|---|---|---|---|---|
| **Grupo Filogenético (Principal)** | {gen_q2['grupo_filogenetico']['auprc']:.4f} | {gen_q2['grupo_filogenetico']['auroc']:.4f} | {gen_q2['grupo_filogenetico']['vme']:.4f} | {gen_q2['grupo_filogenetico']['me']:.4f} | Sem vazamento clonal |
| **Leave-One-Lineage-Out (LOLO)** | {gen_q2['leave_one_lineage_out']['auprc']:.4f} | {gen_q2['leave_one_lineage_out']['auroc']:.4f} | {gen_q2['leave_one_lineage_out']['vme']:.4f} | {gen_q2['leave_one_lineage_out']['me']:.4f} | Linhagens ST inteiras inéditas |
| **Temporal (Treino <= 2022, Teste > 2022)** | {gen_q2['temporal']['auprc']:.4f} | {gen_q2['temporal']['auroc']:.4f} | {gen_q2['temporal']['vme']:.4f} | {gen_q2['temporal']['me']:.4f} | Simula vigilância prospectiva |
| **Geográfica (Países retidos)** | {gen_q2['geografica']['auprc']:.4f} | {gen_q2['geografica']['auroc']:.4f} | {gen_q2['geografica']['vme']:.4f} | {gen_q2['geografica']['me']:.4f} | Robustez a epidemiologia local |
| **Aleatória (Demonstração de Vazamento B5)** | {gen_q2['demonstracao_vazamento_aleatorio_b5']['auprc']:.4f} | {gen_q2['demonstracao_vazamento_aleatorio_b5']['auroc']:.4f} | {gen_q2['demonstracao_vazamento_aleatorio_b5']['vme']:.4f} | {gen_q2['demonstracao_vazamento_aleatorio_b5']['me']:.4f} | Inflação por clones no treino e teste |

#### Confiabilidade, Calibração e Abstenção Conformal (Q4 e REG-11)

- **Brier Score**: `{m_b4['brier_score']:.4f}` | **ECE**: `{m_b4['ece']:.4f}`
- **Cobertura Conformal**: `{m_b4['cobertura_conformal']*100:.1f}%` (Taxa de Abstenção em casos duvidosos: `{m_b4['taxa_abstencao']*100:.1f}%`)
- **Redução de Erros Graves pós-abstenção**: VME cai de `{m_b4['vme']:.4f}` para `{m_b4['vme_pos_abstencao']:.4f}`, e ME cai de `{m_b4['me']:.4f}` para `{m_b4['me_pos_abstencao']:.4f}`.

#### Avaliação de Explicabilidade TreeSHAP Agregada por Locus (Q3, Seção 8.2 e REG-10)

- **Proporção do Top-6 presente no CARD/AMRFinderPlus**: `{xai_q3['proporcao_top_k_no_card']*100:.1f}%` ({', '.join(xai_q3['determinantes_card_recuperados'])})
- **Candidatos novos identificados via k-mers/unitigs (F3)**: `{', '.join(xai_q3['candidatos_novos_top_k']) or 'Nenhum'}`
- **Estabilidade entre dobras (Spearman $\\rho$)**: `{xai_q3['correlacao_estabilidade_entre_dobras']:.4f}`
- **Concordância com GWAS controlado por estrutura populacional**: `{xai_q3['concordancia_gwas_populacional']*100:.1f}%`
- **Teste de remoção in silico do locus principal**: queda média de probabilidade = `-{xai_q3['queda_media_probabilidade_remocao_locus']:.4f}` (Coerente: `{xai_q3['remocao_coerente']}`)

| Locus Genômico | Tipo | Status | Fonte | Importância SHAP Média | Mecanismo Associado |
|---|---|---|---|---|---|
{chr(10).join(linhas_loci)}
"""
        )

    linhas_multiesp = []
    if isinstance(ext6, dict):
        for esp_nome, d_esp in ext6.get("multiespecie", {}).items():
            for farm_e, m_e in d_esp.items():
                linhas_multiesp.append(
                    f"| *{esp_nome}* | `{farm_e}` | `{m_e['versao_modelo']}` | {m_e['auprc_b1']:.4f} | "
                    f"**{m_e['auprc_b4']:.4f}** ({m_e['delta_auprc_vs_b1']:+.4f}) | {m_e['vme_b4']:.4f} | "
                    f"{', '.join(m_e['principais_loci_shap'])} |"
                )

    linhas_val_ext = []
    if isinstance(ext6, dict):
        for farm_v, d_v in ext6.get("validacao_externa_ncbi_pathogen", {}).items():
            linhas_val_ext.append(
                f"| `{farm_v}` | {d_v['auprc_b1_externo']:.4f} | **{d_v['auprc_b4_externo']:.4f}** | "
                f"**{d_v['ganho_auprc_b4_sobre_b1']:+.4f}** | {d_v['vme_b4_externo']:.4f} | {d_v['me_b4_externo']:.4f} |"
            )

    return f"""# Relatório Científico Final — PhyloResistance AI

> **Aviso Ético e Regulatório (REG-15)**: {AVISO_USO_PESQUISA}

---

## 1. Funil de Exclusão e Curadoria de Rótulos (REG-02, REG-06, REG-07)

- **Registros brutos ingeridos**: `{funil['total_registros_brutos']}`
- **Removidos por espécie fora do escopo**: `{funil['removidos_especie_fora_escopo']}`
- **Removidos por fármaco fora do escopo**: `{funil['removidos_farmaco_fora_escopo']}`
- **Removidos por categoria Intermediária (I)**: `{funil['removidos_intermediario']}`
- **Removidos por conflito irreconciliável de rótulo (R vs S)**: `{funil['removidos_conflito_rotulo']}`
- **Removidos por reprovação no QC do genoma (contigs/N50/completude)**: `{funil['removidos_qc_genoma']}`
- **Removidos por desduplicação clonal antes do split (distância Mash <= 0.0005)**: `{funil['removidos_desduplicacao_clonal']}` (`{dedup['total_clones_removidos']}` genoma(s) colapsado(s))
- **Total de registros fenotípicos aprovados para modelagem**: `{funil['total_aprovados_pos_qc_dedup']}`

---

## 2. Resultados do Escopo Principal (*Escherichia coli* — Q1, Q2, Q3 e Q4)

{chr(10).join(secoes_farmacos)}

---

## 3. Extensões da Fase 6: Multi-Espécie (*K. pneumoniae* e *S. aureus*) e Validação Externa

### 3.1 Expansão Multi-Espécie (Seção 3.2 — Fases 2 e 3)

| Espécie | Fármaco | Versão do Modelo | AUPRC Baseline B1 | AUPRC Modelo B4 ($\\Delta$) | VME | Top Loci TreeSHAP |
|---|---|---|---|---|---|---|
{chr(10).join(linhas_multiesp)}

### 3.2 Validação Externa Independente — Padrão NCBI Pathogen Detection AST (Seção 3.1, 7.1 e Fase 6)

| Fármaco (*E. coli*) | AUPRC B1 (Externo) | AUPRC B4 (Externo) | Ganho $\\Delta\\text{{AUPRC}}$ Externo | VME Externo | ME Externo |
|---|---|---|---|---|---|
{chr(10).join(linhas_val_ext)}

---

## 4. Figuras Científicas Geradas (Fase 7)

- `reports/figuras/escada_auprc_ecoli.svg`
- `reports/figuras/calibracao_e_abstencao.svg`
- `reports/figuras/importancia_shap_loci.svg`
- `reports/figuras/extensoes_multiespecie_e_mic.svg`

---

## 5. Limitações e Uso Responsável (Seção 15)

1. **Associação vs. Causalidade (REG-10)**: Os valores TreeSHAP medem contribuição preditiva ao modelo; loci classificados como `candidata` representam hipóteses biológicas para validação experimental.
2. **Domínio de Aplicação**: Os modelos avaliados aplicam-se exclusivamente a isolados das espécies treinadas (*Escherichia coli*, *Klebsiella pneumoniae* e *Staphylococcus aureus*) aprovados pelos critérios de QC da Seção 3.4.
3. **Pesquisa Científica (REG-15)**: Este sistema não é um dispositivo médico e não deve orientar conduta clínica.
"""
