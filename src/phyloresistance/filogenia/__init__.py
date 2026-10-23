"""Módulo de filogenia: distâncias Mash/MinHash, desduplicação clonal (REG-02), clustering e splits anti-vazamento (REG-01)."""

from phyloresistance.filogenia.distancia import (
    RelatorioDesduplicacao,
    calcular_matriz_distancia_mash,
    desduplicar_genomas_quase_identicos,
    inferir_clusters_filogeneticos,
)
from phyloresistance.filogenia.splits import (
    DivisaoTreinoTeste,
    gerar_split_geografico,
    gerar_split_temporal,
    gerar_splits_aleatorios_vazamento,
    gerar_splits_grupo_filogenetico,
    gerar_splits_leave_one_lineage_out,
    verificar_ausencia_vazamento_grupo,
)

__all__ = [
    "DivisaoTreinoTeste",
    "RelatorioDesduplicacao",
    "calcular_matriz_distancia_mash",
    "desduplicar_genomas_quase_identicos",
    "gerar_split_geografico",
    "gerar_split_temporal",
    "gerar_splits_aleatorios_vazamento",
    "gerar_splits_grupo_filogenetico",
    "gerar_splits_leave_one_lineage_out",
    "inferir_clusters_filogeneticos",
    "verificar_ausencia_vazamento_grupo",
]

