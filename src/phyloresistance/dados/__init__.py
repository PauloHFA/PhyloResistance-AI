"""Módulo de ingestão, contratos de dados e curadoria padronizada de rótulos (Seção 3 e REG-06)."""

from phyloresistance.dados.bvbrc import (
    baixar_dados_bvbrc_api,
    calcular_sha256_arquivo,
    carregar_tabela_bvbrc,
    gerar_coorte_especie_extensao,
    gerar_coorte_externa_ncbi_pathogen,
    gerar_coorte_piloto_ecoli,
    gerar_manifesto_dados,
)
from phyloresistance.dados.contratos import (
    FunilExclusao,
    ManifestoVersaoDados,
    RegistroBrutoBVBRC,
    RegistroCurado,
)
from phyloresistance.dados.curadoria import (
    curar_rotulos,
    normalizar_categoria_fenotipo,
)

__all__ = [
    "FunilExclusao",
    "ManifestoVersaoDados",
    "RegistroBrutoBVBRC",
    "RegistroCurado",
    "baixar_dados_bvbrc_api",
    "calcular_sha256_arquivo",
    "carregar_tabela_bvbrc",
    "curar_rotulos",
    "gerar_coorte_especie_extensao",
    "gerar_coorte_externa_ncbi_pathogen",
    "gerar_coorte_piloto_ecoli",
    "gerar_manifesto_dados",
    "normalizar_categoria_fenotipo",
]

