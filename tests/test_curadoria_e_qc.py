"""Testes unitários e de contratos para Curadoria de Rótulos (REG-06), QC (REG-07) e Desduplicação (REG-02)."""

from __future__ import annotations

import numpy as np
import pytest

from phyloresistance.dados.bvbrc import gerar_coorte_piloto_ecoli
from phyloresistance.dados.contratos import RegistroBrutoBVBRC
from phyloresistance.dados.curadoria import curar_rotulos
from phyloresistance.filogenia.distancia import (
    calcular_matriz_distancia_mash,
    desduplicar_genomas_quase_identicos,
)
from phyloresistance.genomica.anotacao import sintetizar_fasta_isolado
from phyloresistance.genomica.qc import (
    ErroValidacaoFASTA,
    avaliar_qc_genoma,
    calcular_estatisticas_fasta,
    validar_fasta_entrada,
)


def test_recalculo_mic_e_tratamento_intermediario_reg06() -> None:
    """Verifica que o MIC tem preferência sobre o rótulo bruto e que 'I' é tratado explicitamente (REG-06)."""
    registros = [
        # Rótulo bruto diz S, mas MIC = 4.0 mg/L para ciprofloxacina (>= 1.0 -> R)
        RegistroBrutoBVBRC(
            genome_id="ISO_01",
            especie="Escherichia coli",
            farmaco="ciprofloxacina",
            fenotipo_bruto="Susceptible",
            mic_mg_l=4.0,
        ),
        # MIC intermediário (0.5 mg/L para ciprofloxacina -> I)
        RegistroBrutoBVBRC(
            genome_id="ISO_02",
            especie="Escherichia coli",
            farmaco="ciprofloxacina",
            fenotipo_bruto="Intermediate",
            mic_mg_l=0.5,
        ),
        # Conflito irreconciliável para ISO_03 (uma leitura R e uma leitura S)
        RegistroBrutoBVBRC(
            genome_id="ISO_03",
            especie="Escherichia coli",
            farmaco="gentamicina",
            fenotipo_bruto="Resistant",
            mic_mg_l=16.0,
        ),
        RegistroBrutoBVBRC(
            genome_id="ISO_03",
            especie="Escherichia coli",
            farmaco="gentamicina",
            fenotipo_bruto="Susceptible",
            mic_mg_l=0.5,
        ),
        # Espécie fora do escopo
        RegistroBrutoBVBRC(
            genome_id="ISO_04",
            especie="Pseudomonas aeruginosa",
            farmaco="ciprofloxacina",
            fenotipo_bruto="Resistant",
            mic_mg_l=8.0,
        ),
    ]

    curados_excluir, funil_excluir = curar_rotulos(registros, modo_intermediario="excluir")
    assert len(curados_excluir) == 1
    assert curados_excluir[0].genome_id == "ISO_01"
    assert curados_excluir[0].fenotipo == "R"
    assert curados_excluir[0].origem_rotulo == "recalculado_mic"
    assert funil_excluir.removidos_intermediario == 1
    assert funil_excluir.removidos_conflito_rotulo == 2
    assert funil_excluir.removidos_especie_fora_escopo == 1

    # Análise de sensibilidade agrupando I com S
    curados_agrupar_s, _ = curar_rotulos(registros, modo_intermediario="agrupar_s")
    assert len(curados_agrupar_s) == 2
    mapa = {c.genome_id: c.fenotipo for c in curados_agrupar_s}
    assert mapa["ISO_02"] == "S"


def test_validacao_fasta_e_qc_genoma_reg07() -> None:
    """Verifica validação de alfabeto ACGTN, cálculo de N50 e gates de QC da Seção 3.4 (REG-07)."""
    fasta_valido = sintetizar_fasta_isolado(
        genome_id="ECOLI_TEST",
        atributos_ativos=["mut_gyrA_S83L", "gene_blaCTX_M_15"],
        n_contigs=40,
        tamanho_total_pb=4_900_000,
    )
    stats = calcular_estatisticas_fasta(fasta_valido)
    assert stats["contigs"] == 40
    assert 4_500_000 <= stats["tamanho_pb"] <= 5_500_000
    assert stats["n50"] >= 50_000

    qc_ok = avaliar_qc_genoma(
        genome_id="ECOLI_TEST",
        tamanho_pb=int(stats["tamanho_pb"]),
        contigs=int(stats["contigs"]),
        n50=int(stats["n50"]),
        fracao_n=float(stats["fracao_n"]),
    )
    assert qc_ok.aprovado is True
    assert len(qc_ok.motivos_reprovacao) == 0

    # Genoma fragmentado ou de espécie divergente deve reprovar no QC
    qc_falho = avaliar_qc_genoma(
        genome_id="ECOLI_BAD",
        tamanho_pb=2_800_000,
        contigs=950,
        n50=5_000,
        especie_declarada="Escherichia coli",
        especie_inferida="Klebsiella pneumoniae",
    )
    assert qc_falho.aprovado is False
    assert len(qc_falho.motivos_reprovacao) >= 3

    # Sequência com caracteres inválidos fora de ACGTN deve levantar ErroValidacaoFASTA
    with pytest.raises(ErroValidacaoFASTA):
        validar_fasta_entrada(">contig_1\nACGTXYZ123\n")


def test_desduplicacao_clones_quase_identicos_reg02() -> None:
    """Verifica que clones idênticos/quase-idênticos são detectados e removidos antes do split (REG-02)."""
    X = np.array(
        [
            [1, 1, 0, 1, 0, 1],
            [1, 1, 0, 1, 0, 1],  # clone exato do isolado 0
            [0, 0, 1, 0, 1, 0],
            [1, 0, 1, 1, 0, 0],
        ],
        dtype=np.float64,
    )
    ids = ["ISO_A", "ISO_A_CLONE", "ISO_B", "ISO_C"]
    dist = calcular_matriz_distancia_mash(X, k=31)
    assert dist[0, 1] == 0.0

    idx_mantidos, relatorio = desduplicar_genomas_quase_identicos(
        genome_ids=ids, matriz_distancia=dist, limiar_distancia=0.0005
    )
    assert idx_mantidos == [0, 2, 3]
    assert relatorio.total_clones_removidos == 1
    assert "ISO_A_CLONE" in relatorio.ids_removidos


def test_coorte_piloto_exercita_funil_completo() -> None:
    """Verifica que a coorte piloto contém casos reais de exclusão e mantém isolados balanceados."""
    registros, df_meta, df_geno = gerar_coorte_piloto_ecoli(n_isolados_base=120, seed=42)
    curados, funil = curar_rotulos(registros, modo_intermediario="excluir")
    assert funil.removidos_especie_fora_escopo >= 1
    assert funil.removidos_farmaco_fora_escopo >= 1
    assert funil.removidos_intermediario >= 1
    assert funil.removidos_conflito_rotulo >= 2
    assert len(curados) >= 300

