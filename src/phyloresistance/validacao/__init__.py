"""Módulo de validação: métricas clínicas/microbiológicas (VME, ME, AUPRC, Brier) e bootstrap em blocos por cluster."""

from phyloresistance.validacao.metricas import (
    ComparacaoPareadaBaseline,
    IntervaloConfiancaBloco,
    ResultadoMetricasModelo,
    bootstrap_em_blocos_por_cluster,
    calcular_ece,
    calcular_metricas_completas,
    comparar_modelo_vs_baseline_pareado,
)

__all__ = [
    "ComparacaoPareadaBaseline",
    "IntervaloConfiancaBloco",
    "ResultadoMetricasModelo",
    "bootstrap_em_blocos_por_cluster",
    "calcular_ece",
    "calcular_metricas_completas",
    "comparar_modelo_vs_baseline_pareado",
]

