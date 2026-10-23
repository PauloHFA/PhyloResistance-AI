"""Módulo de genômica: controle de qualidade (QC), k-mers canônicos, unitigs e anotação AMR (F1, F2, F3)."""

from phyloresistance.genomica.anotacao import (
    anotar_isolado_fasta,
    construir_matrizes_atributos_coorte,
    sintetizar_fasta_isolado,
)
from phyloresistance.genomica.fastq import (
    ErroValidacaoFASTQ,
    ResultadoTrimagemMontagemFASTQ,
    ler_registros_fastq,
    sintetizar_fastq_isolado,
    trimar_e_montar_fastq_para_fasta,
)
from phyloresistance.genomica.kmers import (
    FiltroSelecaoKmersDobra,
    agrupar_kmers_equivalentes_unitigs,
    complemento_reverso,
    extrair_frequencias_kmers_curtos,
    extrair_kmers_canonicos,
)
from phyloresistance.genomica.qc import (
    ErroValidacaoFASTA,
    ResultadoQCGenoma,
    avaliar_qc_coorte,
    avaliar_qc_genoma,
    calcular_estatisticas_fasta,
    validar_fasta_entrada,
)

__all__ = [
    "ErroValidacaoFASTA",
    "ErroValidacaoFASTQ",
    "FiltroSelecaoKmersDobra",
    "ResultadoQCGenoma",
    "ResultadoTrimagemMontagemFASTQ",
    "agrupar_kmers_equivalentes_unitigs",
    "anotar_isolado_fasta",
    "avaliar_qc_coorte",
    "avaliar_qc_genoma",
    "calcular_estatisticas_fasta",
    "complemento_reverso",
    "construir_matrizes_atributos_coorte",
    "extrair_frequencias_kmers_curtos",
    "extrair_kmers_canonicos",
    "ler_registros_fastq",
    "sintetizar_fasta_isolado",
    "sintetizar_fastq_isolado",
    "trimar_e_montar_fastq_para_fasta",
    "validar_fasta_entrada",
]

