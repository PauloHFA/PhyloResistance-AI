"""Módulo de modelagem: escada de modelos B0 a B5, controle k-curto, Optuna aninhado, calibração e abstenção conformal."""

from phyloresistance.modelos.calibracao import (
    CalibradorConformalAMR,
    PredicaoComIncerteza,
)
from phyloresistance.modelos.escada import (
    BaselineRegrasB1,
    ClassificadorMajoritarioB0,
    ModeloArvoreCuradoB3,
    ModeloControleKmersCurtos,
    ModeloLinearEsparsoKmersB5,
    ModeloLogisticoRegularizadoB2,
    PipelineCompletoB4,
    calcular_peso_classe_balanceado,
    otimizar_xgb_optuna_por_grupo,
)

__all__ = [
    "BaselineRegrasB1",
    "CalibradorConformalAMR",
    "ClassificadorMajoritarioB0",
    "ModeloArvoreCuradoB3",
    "ModeloControleKmersCurtos",
    "ModeloLinearEsparsoKmersB5",
    "ModeloLogisticoRegularizadoB2",
    "PipelineCompletoB4",
    "PredicaoComIncerteza",
    "calcular_peso_classe_balanceado",
    "otimizar_xgb_optuna_por_grupo",
]

