"""Orquestrador ponta a ponta do pipeline PhyloResistance AI (Etapas 1 a 8) e Registro de Modelos."""

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
    FARMACOS_PILOTO,
    GLOBAL_RANDOM_SEED,
    VERSOES_BANCOS,
    VERSOES_MODELOS,
)
from phyloresistance.dados.bvbrc import (
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
from phyloresistance.validacao.metricas import (
    ComparacaoPareadaBaseline,
    ResultadoMetricasModelo,
    calcular_metricas_completas,
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


class RegistroModelosServico:
    """Gerencia os modelos treinados por espécie e fármaco para o serviço FastAPI."""

    def __init__(self) -> None:
        self.modelos_por_farmaco: dict[str, ArtefatoModeloFarmaco] = {}
        self.especies_suportadas: set[str] = {ESPECIE_PILOTO}
        self.versoes_bancos: dict[str, str] = dict(VERSOES_BANCOS)

    def esta_pronto(self) -> bool:
        return len(self.modelos_por_farmaco) > 0

    def registrar(self, artefato: ArtefatoModeloFarmaco) -> None:
        self.modelos_por_farmaco[artefato.farmaco] = artefato

    def salvar(self, caminho: str | Path) -> None:
        p = Path(caminho)
        p.parent.mkdir(parents=True, exist_ok=True)
        with open(p, "wb") as f:
            pickle.dump(self.modelos_por_farmaco, f)

    def carregar(self, caminho: str | Path) -> None:
        p = Path(caminho)
        with open(p, "rb") as f:
            self.modelos_por_farmaco = pickle.load(f)

    def predizer_isolado_fasta(
        self,
        conteudo_ou_caminho_fasta: str | Path,
        farmacos_solicitados: list[str] | None = None,
    ) -> list[dict[str, object]]:
        """
        Realiza a inferência completa para um isolado FASTA:
        1) Anotação genômica F1 + F2 + F3
        2) Predição com modelo B4 + calibração conformal e abstenção
        3) Explicação TreeSHAP agregada por locus e cruzada com o CARD/AMRFinderPlus
        """
        atributos_detectados = anotar_isolado_fasta(conteudo_ou_caminho_fasta)
        farmacos = farmacos_solicitados or list(self.modelos_por_farmaco.keys())
        predicoes_saida: list[dict[str, object]] = []

        for farmaco in farmacos:
            if farmaco not in self.modelos_por_farmaco:
                predicoes_saida.append(
                    {
                        "farmaco": farmaco,
                        "status": "nao_suportado",
                        "mensagem": f"Fármaco '{farmaco}' fora do domínio treinado ({list(self.modelos_por_farmaco.keys())}).",
                    }
                )
                continue

            art = self.modelos_por_farmaco[farmaco]
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


def executar_pipeline_completo(
    output_dir: str | Path = "data/processed",
    reports_dir: str | Path = "reports",
    n_isolados_base: int = 360,
    n_boot: int = 200,
    seed: int = GLOBAL_RANDOM_SEED,
) -> tuple[RegistroModelosServico, dict[str, object]]:
    """
    Executa de ponta a ponta as Etapas 1 a 8 da arquitetura do PhyloResistance AI:
    [1] Ingestão e curadoria padronizada de rótulos (mais análise de sensibilidade de 'I')
    [2] QC dos genomas + confirmação de espécie + desduplicação clonal (REG-02, REG-07)
    [3] Estrutura filogenética (distâncias Mash -> clusters / ST)
    [4] Engenharia de atributos: F1 (genes AMR), F2 (mutações pontuais), F3 (k-mers k=31 esparsos) e Controle k=6
    [5] Divisão por grupo filogenético (principal), LOLO, Temporal, Geográfica e Aleatória (demonstração de vazamento)
    [6] Escada de modelos (B0, B1, B2, B3, B4, B5 e Controle k=6)
    [7] Calibração probabilística + incerteza / abstenção conformal (Q4, REG-11)
    [8] Explicabilidade TreeSHAP agregada por locus + 4 testes de avaliação da Seção 8.2 (Q3, REG-10)
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
    # Análise de sensibilidade (Seção 3.3: agrupar I com S)
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

    # Filtra matrizes e tabelas mantendo APENAS isolados aprovados no QC e únicos pós-desduplicação
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

    # Salvar artefatos das Fases 1 e 2 em data/processed/
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

    # --- ETAPAS 5, 6, 7 e 8: Modelagem, Validação, Calibração e XAI por Fármaco ---
    nomes_f1_f2 = [c for c in df_f1_f2.columns if c != "genome_id"]
    X_f1_f2_full = df_f1_f2[nomes_f1_f2].to_numpy(dtype=np.float64)
    mapa_id_para_idx = {gid: idx for idx, gid in enumerate(df_meta_unico["genome_id"].values)}

    registro_servico = RegistroModelosServico()
    resultados_por_farmaco: dict[str, object] = {}

    for farmaco in FARMACOS_PILOTO:
        df_sub = df_curados_final[df_curados_final["farmaco"] == farmaco].copy()
        # Alinha na mesma ordem de df_meta_unico
        indices_alinhados = [mapa_id_para_idx[gid] for gid in df_sub["genome_id"].values]
        y = df_sub["y"].to_numpy(dtype=np.int32)
        ids_isolados = list(df_sub["genome_id"].values)
        grupos = np.array([clusters_filogeneticos[i] for i in indices_alinhados], dtype=str)
        anos = df_meta_unico["ano"].to_numpy(dtype=np.int32)[indices_alinhados]
        paises = df_meta_unico["pais"].to_numpy(dtype=str)[indices_alinhados]

        X_f12 = X_f1_f2_full[indices_alinhados]
        X_f3 = X_f3_sparse[indices_alinhados]
        X_k6 = X_k6_controle[indices_alinhados]

        # 5.1) Validação Principal por Grupo Filogenético (StratifiedGroupKFold, 5 dobras — REG-01)
        splits_grupo = gerar_splits_grupo_filogenetico(y, grupos, n_splits=5, seed=seed)
        n = len(y)
        oof_probs: dict[str, np.ndarray] = {
            "B0_Classe_Majoritaria": np.zeros(n, dtype=np.float64),
            "B1_Regras_CARD_AMRFinder": np.zeros(n, dtype=np.float64),
            "B2_LogReg_L1_F1_F2": np.zeros(n, dtype=np.float64),
            "B3_XGBoost_F1_F2": np.zeros(n, dtype=np.float64),
            "B4_XGBoost_F1_F2_F3": np.zeros(n, dtype=np.float64),
            "B5_Linear_Esparso_F3": np.zeros(n, dtype=np.float64),
            "Controle_Kmers_Curtos_k6": np.zeros(n, dtype=np.float64),
        }
        oof_probs_b4_brutas = np.zeros(n, dtype=np.float64)
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

            # Coleta ranking SHAP agregado por locus nesta dobra para medir Estabilidade (Seção 8.2)
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

            m_ctrl = ModeloControleKmersCurtos(seed=seed).fit(X_k6[tr], y[tr])
            oof_probs["Controle_Kmers_Curtos_k6"][te] = m_ctrl.predict_proba(X_k6[te])[:, 1]

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
                n_boot=n_boot,
                seed=seed,
            )
            metricas_escada[nome_mod] = res_m.model_dump()

        # Comparação pareada com IC 95% por cluster contra B1 (Q1 e REG-04)
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

        # 5.2) Validações de Generalização (Q2 e REG-09): LOLO, Temporal, Geográfica e Split Aleatório (Vazamento)
        # a) Split Temporal
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

        # b) Split Geográfico
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

        # c) Leave-One-Lineage-Out (LOLO)
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

        # d) Demonstração do Vazamento Filogenético no Split Aleatório (Seção 7.1)
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

        # Treina modelo final B4 no conjunto completo para registro no serviço e avaliação XAI
        modelo_b4_final = PipelineCompletoB4(
            nomes_f1_f2=nomes_f1_f2,
            nomes_f3=nomes_f3,
            max_kmers_f3=50,
            seed=seed,
        ).fit(X_f12, X_f3, y, ids_treino=ids_isolados)

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
        }

    # Salvar registro de modelos para a API
    caminho_modelos = dir_proc / "modelos_servico.pkl"
    registro_servico.salvar(caminho_modelos)

    # Gerar Manifesto de Hashes SHA-256 (Fase 0 / REG-13)
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
        "farmacos_piloto": list(FARMACOS_PILOTO),
        "versoes_bancos": dict(VERSOES_BANCOS),
        "funil_exclusao": funil.model_dump(),
        "analise_sensibilidade_intermediario": {
            "total_excluindo_I": len(curados_principal),
            "total_agrupando_I_com_S": len(curados_sens_i_como_s),
        },
        "desduplicacao": relatorio_dedup.model_dump(),
        "resultados_por_farmaco": resultados_por_farmaco,
        "manifesto_dados": manifesto.model_dump(),
        "aviso": AVISO_USO_PESQUISA,
    }

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
    """Gera o relatório científico final (formato Quarto / Markdown) respondendo Q1, Q2, Q3 e Q4."""
    funil = rel["funil_exclusao"]  # type: ignore[index]
    dedup = rel["desduplicacao"]  # type: ignore[index]
    res_farmacos = rel["resultados_por_farmaco"]  # type: ignore[index]

    secoes_farmacos: list[str] = []
    for farmaco, dados in res_farmacos.items():  # type: ignore[union-attr]
        escada = dados["escada_modelos_grupo_filogenetico"]
        comp_b4 = dados["comparacao_pareada_vs_b1"]["B4_vs_B1"]
        gen_q2 = dados["generalizacao_q2"]
        xai_q3 = dados["explicabilidade_q3"]

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

## 2. Resultados por Fármaco (Perguntas Científicas Q1, Q2, Q3 e Q4)

{chr(10).join(secoes_farmacos)}

---

## 3. Limitações e Uso Responsável (Seção 15)

1. **Associação vs. Causalidade (REG-10)**: Os valores TreeSHAP medem contribuição preditiva ao modelo; loci classificados como `candidata` representam hipóteses biológicas para validação experimental (ex.: mutagênese sítio-dirigida ou ensaio funcional).
2. **Domínio de Aplicação**: Os modelos avaliados aplicam-se exclusivamente a isolados montados de *Escherichia coli* aprovados pelos critérios de QC da Seção 3.4 para os fármacos `ciprofloxacina`, `gentamicina` e `cefuroxima`.
3. **Pesquisa Científica (REG-15)**: Este sistema não é um dispositivo médico e não deve orientar conduta clínica.
"""
