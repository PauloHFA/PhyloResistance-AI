"""Módulo de validação: métricas clínicas/microbiológicas (VME, ME, AUPRC, Brier) e bootstrap em blocos por cluster."""

from phyloresistance.validacao.metricas import (
    ComparacaoPareadaBaseline,
    IntervaloConfiancaBloco,
    ResultadoAvaliacaoMIC,
    ResultadoMetricasModelo,
    bootstrap_em_blocos_por_cluster,
    calcular_ece,
    calcular_metricas_completas,
    calcular_metricas_mic_iso20776,
    comparar_modelo_vs_baseline_pareado,
)

__all__ = [
    "ComparacaoPareadaBaseline",
    "IntervaloConfiancaBloco",
    "ResultadoAvaliacaoMIC",
    "ResultadoMetricasModelo",
    "bootstrap_em_blocos_por_cluster",
    "calcular_ece",
    "calcular_metricas_completas",
    "calcular_metricas_mic_iso20776",
    "comparar_modelo_vs_baseline_pareado",
]

