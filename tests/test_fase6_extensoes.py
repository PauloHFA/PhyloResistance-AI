"""Testes automatizados das extensões da Fase 6: FASTQ, Regressão log2(MIC) ISO 20776-2, B6, Multi-Saída, Multi-Espécie e Validação Externa."""

from __future__ import annotations

import numpy as np
import pytest

from phyloresistance.config import PONTOS_CORTE_ECOLI, PONTOS_CORTE_SAUREUS
from phyloresistance.dados.bvbrc import (
    gerar_coorte_especie_extensao,
    gerar_coorte_externa_ncbi_pathogen,
    gerar_coorte_piloto_ecoli,
)
from phyloresistance.dados.curadoria import (
    converter_curados_para_dataframe,
    curar_rotulos,
)
from phyloresistance.genomica.anotacao import (
    anotar_isolado_fasta,
    construir_matrizes_atributos_coorte,
)
from phyloresistance.genomica.fastq import (
    ErroValidacaoFASTQ,
    sintetizar_fastq_isolado,
    trimar_e_montar_fastq_para_fasta,
)
from phyloresistance.modelos.escada import BaselineRegrasB1, PipelineCompletoB4
from phyloresistance.modelos.extensoes import (
    ModeloMultiSaidaAMR,
    ModeloSequenciaConvolucionalB6,
    RegressorMICLog2,
)
from phyloresistance.validacao.metricas import (
    calcular_metricas_completas,
    calcular_metricas_mic_iso20776,
)


def test_pipeline_fastq_trimagem_e_montagem() -> None:
    """Verifica trimagem de leituras Phred e montagem FASTQ -> FASTA preservando determinantes AMR."""
    fastq_str = sintetizar_fastq_isolado(
        genome_id="FASTQ_ISO_01",
        atributos_ativos=["mut_gyrA_S83L", "gene_sa_mecA"],
        n_reads_alta_qualidade=35,
        n_reads_baixa_qualidade=6,
    )
    res_montagem = trimar_e_montar_fastq_para_fasta(
        fastq_str,
        genome_id="FASTQ_ISO_01",
        min_phred=20.0,
        tamanho_alvo_montagem_pb=2_900_000,
    )
    assert res_montagem.total_reads_brutos == 41
    # As 6 leituras de baixa qualidade (Phred 8) devem ter sido descartadas na trimagem
    assert res_montagem.total_reads_pos_trimagem == 35
    assert res_montagem.phred_medio >= 35.0
    assert res_montagem.contigs_montados == 40
    assert res_montagem.n50_montado >= 25_000

    atributos = anotar_isolado_fasta(res_montagem.fasta_montado)
    assert atributos["mut_gyrA_S83L"] == 1
    assert atributos["gene_sa_mecA"] == 1

    with pytest.raises(ErroValidacaoFASTQ):
        trimar_e_montar_fastq_para_fasta("linha_invalida_sem_arroba\nACGT\n+\nIIII\n")


def test_regressao_mic_log2_e_concordancia_essencial_iso20776() -> None:
    """Verifica regressão de log2(MIC) e Concordância Essencial (+-1 diluição dupla, ISO 20776-2)."""
    regs, _, df_geno = gerar_coorte_piloto_ecoli(n_isolados_base=140, seed=42)
    df_geno_base = df_geno.iloc[:140].reset_index(drop=True)
    df_f12, X_f3, _, _ = construir_matrizes_atributos_coorte(df_geno_base, seed=42)
    nomes_f12 = [c for c in df_f12.columns if c != "genome_id"]
    X_f12 = df_f12[nomes_f12].to_numpy(dtype=np.float64)

    mapa_mic = {
        r.genome_id: float(r.mic_mg_l)
        for r in regs
        if r.farmaco == "ciprofloxacina" and r.genome_id.startswith("ECOLI_1") and r.mic_mg_l is not None
    }
    mics = np.array([mapa_mic[gid] for gid in df_geno_base["genome_id"].values], dtype=np.float64)

    reg_mic = RegressorMICLog2(max_kmers_f3=40, seed=42).fit(X_f12[:100], X_f3[:100], mics[:100])
    log2_pred = reg_mic.predict_log2_mic(X_f12[100:], X_f3[100:])
    mics_pred = reg_mic.predict_mic_mg_l(X_f12[100:], X_f3[100:])

    assert len(mics_pred) == 40
    assert np.all(mics_pred > 0)

    aval = calcular_metricas_mic_iso20776(
        farmaco="ciprofloxacina",
        mic_real_mg_l=mics[100:],
        log2_mic_pred=log2_pred,
        ponto_corte=PONTOS_CORTE_ECOLI["ciprofloxacina"],
    )
    assert aval.concordancia_essencial_1_diluicao >= 0.70
    assert aval.concordancia_categorica >= 0.75


def test_modelo_sequencia_b6_e_modelo_multisaida() -> None:
    """Verifica o modelo convolucional de sequência nível B6 e o classificador multi-saída."""
    regs, _, df_geno = gerar_coorte_piloto_ecoli(n_isolados_base=120, seed=42)
    df_geno_base = df_geno.iloc[:120].reset_index(drop=True)
    df_f12, X_f3, _, _ = construir_matrizes_atributos_coorte(df_geno_base, seed=42)
    nomes_f12 = [c for c in df_f12.columns if c != "genome_id"]
    X_f12 = df_f12[nomes_f12].to_numpy(dtype=np.float64)

    curados, _ = curar_rotulos(regs)
    df_cur = converter_curados_para_dataframe(curados)
    pivo = df_cur[df_cur["genome_id"].isin(df_geno_base["genome_id"])].pivot(
        index="genome_id", columns="farmaco", values="y"
    ).loc[df_geno_base["genome_id"].values]
    Y_multi = pivo[["ciprofloxacina", "gentamicina", "cefuroxima"]].to_numpy(dtype=np.int32)

    b6 = ModeloSequenciaConvolucionalB6(seed=42).fit(X_f12[:90], X_f3[:90], Y_multi[:90, 0])
    probs_b6 = b6.predict_proba(X_f12[90:], X_f3[90:])[:, 1]
    assert probs_b6.shape == (30,)
    assert np.all((probs_b6 >= 0.0) & (probs_b6 <= 1.0))

    multi = ModeloMultiSaidaAMR(seed=42).fit(X_f12[:90], Y_multi[:90])
    probs_mat = multi.predict_proba_matriz(X_f12[90:])
    assert probs_mat.shape == (30, 3)


def test_multiespecie_kpneumoniae_e_saureus_mrsa() -> None:
    """Verifica geração, curadoria e modelagem para K. pneumoniae (Fase 2) e S. aureus MRSA (Fase 3)."""
    for especie, farmaco_teste in [
        ("Klebsiella pneumoniae", "ceftriaxona"),
        ("Staphylococcus aureus", "oxacilina"),
    ]:
        regs, _, df_geno = gerar_coorte_especie_extensao(especie=especie, n_isolados_base=100, seed=42)
        curados, _ = curar_rotulos(regs, especie_alvo=especie, farmacos_alvo=(farmaco_teste,))
        assert len(curados) == 100
        if especie == "Staphylococcus aureus":
            assert PONTOS_CORTE_SAUREUS["oxacilina"].classificar(8.0) == "R"
            assert PONTOS_CORTE_SAUREUS["oxacilina"].classificar(0.5) == "S"
        df_f12, X_f3, nomes_f3, _ = construir_matrizes_atributos_coorte(df_geno, especie=especie, seed=42)
        assert len(df_f12) == 100
        assert X_f3.shape[0] == 100
        assert len(nomes_f3) > 10


def test_validacao_externa_ncbi_pathogen_detection() -> None:
    """
    Critério de Aceite da Fase 6:
    Verifica ganho medido do modelo B4 sobre o Baseline B1 em coorte de validação externa independente (NCBI Pathogen Detection).
    """
    regs_tr, _, df_geno_tr = gerar_coorte_piloto_ecoli(n_isolados_base=160, seed=42)
    df_geno_tr = df_geno_tr.iloc[:160].reset_index(drop=True)
    df_f12_tr, X_f3_tr, nomes_f3, _ = construir_matrizes_atributos_coorte(df_geno_tr, seed=42)
    nomes_f12 = [c for c in df_f12_tr.columns if c != "genome_id"]
    X_f12_tr_arr = df_f12_tr[nomes_f12].to_numpy(dtype=np.float64)
    mapa_y_tr = {
        r.genome_id: (1 if r.fenotipo_bruto == "Resistant" else 0)
        for r in regs_tr
        if r.farmaco == "ciprofloxacina" and r.genome_id.startswith("ECOLI_1")
    }
    y_tr = np.array([mapa_y_tr[gid] for gid in df_geno_tr["genome_id"].values], dtype=np.int32)

    b1 = BaselineRegrasB1("ciprofloxacina", nomes_f12).fit(X_f12_tr_arr, y_tr)
    b4 = PipelineCompletoB4(nomes_f12, nomes_f3, max_kmers_f3=45, seed=42).fit(
        X_f12_tr_arr, X_f3_tr, y_tr
    )

    # Coorte externa independente NCBI Pathogen Detection
    regs_ext, _, df_geno_ext = gerar_coorte_externa_ncbi_pathogen(
        especie="Escherichia coli", n_isolados=100, seed=2099
    )
    df_f12_ext, X_f3_ext, _, _ = construir_matrizes_atributos_coorte(df_geno_ext, seed=2099)
    X_f12_ext_arr = df_f12_ext[nomes_f12].to_numpy(dtype=np.float64)
    mapa_y_ext = {
        r.genome_id: (1 if r.fenotipo_bruto == "Resistant" else 0)
        for r in regs_ext
        if r.farmaco == "ciprofloxacina"
    }
    y_ext = np.array([mapa_y_ext[gid] for gid in df_geno_ext["genome_id"].values], dtype=np.int32)

    met_b1_ext = calcular_metricas_completas(
        "B1", "externa", y_ext, b1.predict_proba(X_f12_ext_arr)[:, 1], n_boot=0
    )
    met_b4_ext = calcular_metricas_completas(
        "B4", "externa", y_ext, b4.predict_proba(X_f12_ext_arr, X_f3_ext)[:, 1], n_boot=0
    )

    assert met_b4_ext.auprc > met_b1_ext.auprc
    assert met_b4_ext.vme < met_b1_ext.vme

