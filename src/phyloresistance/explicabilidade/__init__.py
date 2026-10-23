"""Módulo de explicabilidade (XAI): TreeSHAP, agregação por locus, cruzamento CARD e testes da Seção 8.2."""

from phyloresistance.explicabilidade.mapeamento_genomico import (
    AlinhamentoKmerLocus,
    LocusReferenciaGFF,
    alinhar_kmers_contra_referencia_gff3,
    carregar_anotacao_gff3,
    gerar_referencia_e_gff3_sintetico_ecoli,
)
from phyloresistance.explicabilidade.shap_locus import (
    AvaliacaoExplicabilidadeQ3,
    EvidenciaGenomica,
    InfoLocusGenomico,
    agregar_shap_por_locus,
    avaliar_concordancia_gwas_populacional,
    avaliar_estabilidade_explicacoes,
    avaliar_explicabilidade_completa,
    avaliar_recuperacao_determinantes_card,
    avaliar_remocao_locus_principal,
    calcular_valores_treeshap,
    construir_mapa_atributo_para_locus,
    extrair_evidencias_isolado,
    gerar_relatorio_html_isolado,
    gerar_relatorio_pdf_isolado,
)

__all__ = [
    "AlinhamentoKmerLocus",
    "AvaliacaoExplicabilidadeQ3",
    "EvidenciaGenomica",
    "InfoLocusGenomico",
    "LocusReferenciaGFF",
    "agregar_shap_por_locus",
    "alinhar_kmers_contra_referencia_gff3",
    "avaliar_concordancia_gwas_populacional",
    "avaliar_estabilidade_explicacoes",
    "avaliar_explicabilidade_completa",
    "avaliar_recuperacao_determinantes_card",
    "avaliar_remocao_locus_principal",
    "calcular_valores_treeshap",
    "carregar_anotacao_gff3",
    "construir_mapa_atributo_para_locus",
    "extrair_evidencias_isolado",
    "gerar_referencia_e_gff3_sintetico_ecoli",
    "gerar_relatorio_html_isolado",
    "gerar_relatorio_pdf_isolado",
]

