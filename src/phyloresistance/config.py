"""Configurações globais, sementes, pontos de corte clínicos e catálogos CARD/AMRFinderPlus multi-espécie."""

from __future__ import annotations

from dataclasses import dataclass

GLOBAL_RANDOM_SEED: int = 42
AVISO_USO_PESQUISA: str = "Uso em pesquisa. Não substitui o antibiograma."

VERSOES_BANCOS: dict[str, str] = {
    "amrfinderplus": "3.12.8",
    "card": "3.2.9",
    "refseq_ecoli": "NC_000913.3",
    "refseq_kpneumoniae": "NC_016845.1",
    "refseq_saureus": "NC_007795.1",
}

ESPECIE_PILOTO: str = "Escherichia coli"
FARMACOS_PILOTO: tuple[str, ...] = ("ciprofloxacina", "gentamicina", "cefuroxima")

ESPECIES_SUPORTADAS: tuple[str, ...] = (
    "Escherichia coli",
    "Klebsiella pneumoniae",
    "Staphylococcus aureus",
)

FARMACOS_POR_ESPECIE: dict[str, tuple[str, ...]] = {
    "Escherichia coli": ("ciprofloxacina", "gentamicina", "cefuroxima"),
    "Klebsiella pneumoniae": ("ciprofloxacina", "gentamicina", "ceftriaxona"),
    "Staphylococcus aureus": ("oxacilina", "ciprofloxacina"),
}

VERSOES_MODELOS: dict[str, str] = {
    "ciprofloxacina": "ecoli-cip-1.2.0",
    "gentamicina": "ecoli-gen-1.2.0",
    "cefuroxima": "ecoli-cxm-1.2.0",
}

VERSOES_MODELOS_POR_ESPECIE: dict[str, dict[str, str]] = {
    "Escherichia coli": {
        "ciprofloxacina": "ecoli-cip-1.2.0",
        "gentamicina": "ecoli-gen-1.2.0",
        "cefuroxima": "ecoli-cxm-1.2.0",
    },
    "Klebsiella pneumoniae": {
        "ciprofloxacina": "kpneu-cip-1.2.0",
        "gentamicina": "kpneu-gen-1.2.0",
        "ceftriaxona": "kpneu-cro-1.2.0",
    },
    "Staphylococcus aureus": {
        "oxacilina": "saureus-oxa-1.2.0",
        "ciprofloxacina": "saureus-cip-1.2.0",
    },
}

# Colunas de metadados epidemiológicos que NUNCA podem entrar como atributos no modelo (REG-08)
COLUNAS_METADADOS_PROIBIDAS: frozenset[str] = frozenset(
    {
        "genome_id",
        "especie",
        "pais",
        "ano",
        "st",
        "linhagem",
        "cluster_filogenetico",
        "metodo_laboratorio",
        "padrao",
        "versao_ponto_corte",
        "mic_mg_l",
        "fenotipo",
        "y",
    }
)


@dataclass(frozen=True)
class PontoCorteMIC:
    """Pontos de corte clínicos CLSI/EUCAST em mg/L para recálculo determinístico de S/I/R."""

    farmaco: str
    suscetivel_max: float
    resistente_min: float
    padrao: str = "CLSI M100 / EUCAST"
    versao: str = "2024.1"

    def classificar(self, mic_mg_l: float) -> str:
        if mic_mg_l <= self.suscetivel_max:
            return "S"
        if mic_mg_l >= self.resistente_min:
            return "R"
        return "I"


PONTOS_CORTE_ECOLI: dict[str, PontoCorteMIC] = {
    "ciprofloxacina": PontoCorteMIC(
        farmaco="ciprofloxacina",
        suscetivel_max=0.25,
        resistente_min=1.0,
        padrao="CLSI/EUCAST",
        versao="2024.1",
    ),
    "gentamicina": PontoCorteMIC(
        farmaco="gentamicina",
        suscetivel_max=2.0,
        resistente_min=8.0,
        padrao="CLSI/EUCAST",
        versao="2024.1",
    ),
    "cefuroxima": PontoCorteMIC(
        farmaco="cefuroxima",
        suscetivel_max=8.0,
        resistente_min=32.0,
        padrao="CLSI",
        versao="2024.1",
    ),
    "ceftriaxona": PontoCorteMIC(
        farmaco="ceftriaxona",
        suscetivel_max=1.0,
        resistente_min=4.0,
        padrao="CLSI/EUCAST",
        versao="2024.1",
    ),
}

PONTOS_CORTE_KPNEUMONIAE: dict[str, PontoCorteMIC] = {
    "ciprofloxacina": PontoCorteMIC(
        farmaco="ciprofloxacina",
        suscetivel_max=0.25,
        resistente_min=1.0,
        padrao="CLSI/EUCAST",
        versao="2024.1",
    ),
    "gentamicina": PontoCorteMIC(
        farmaco="gentamicina",
        suscetivel_max=2.0,
        resistente_min=8.0,
        padrao="CLSI/EUCAST",
        versao="2024.1",
    ),
    "ceftriaxona": PontoCorteMIC(
        farmaco="ceftriaxona",
        suscetivel_max=1.0,
        resistente_min=4.0,
        padrao="CLSI/EUCAST",
        versao="2024.1",
    ),
}

PONTOS_CORTE_SAUREUS: dict[str, PontoCorteMIC] = {
    "oxacilina": PontoCorteMIC(
        farmaco="oxacilina",
        suscetivel_max=2.0,
        resistente_min=4.0,
        padrao="CLSI M100 (MRSA)",
        versao="2024.1",
    ),
    "ciprofloxacina": PontoCorteMIC(
        farmaco="ciprofloxacina",
        suscetivel_max=1.0,
        resistente_min=4.0,
        padrao="CLSI M100",
        versao="2024.1",
    ),
}

PONTOS_CORTE_POR_ESPECIE: dict[str, dict[str, PontoCorteMIC]] = {
    "Escherichia coli": PONTOS_CORTE_ECOLI,
    "Klebsiella pneumoniae": PONTOS_CORTE_KPNEUMONIAE,
    "Staphylococcus aureus": PONTOS_CORTE_SAUREUS,
}


@dataclass(frozen=True)
class CriterioQCEspecie:
    """Limiares de controle de qualidade da montagem por espécie (Seção 3.4 e REG-07)."""

    especie: str
    tamanho_min_pb: int = 4_400_000
    tamanho_max_pb: int = 5_800_000
    max_contigs: int = 500
    min_n50_pb: int = 20_000
    max_fracao_n: float = 0.02
    min_completude: float = 95.0
    max_contaminacao: float = 5.0


CRITERIOS_QC_POR_ESPECIE: dict[str, CriterioQCEspecie] = {
    "Escherichia coli": CriterioQCEspecie(
        especie="Escherichia coli",
        tamanho_min_pb=4_400_000,
        tamanho_max_pb=5_800_000,
        max_contigs=500,
        min_n50_pb=20_000,
        max_fracao_n=0.02,
        min_completude=95.0,
        max_contaminacao=5.0,
    ),
    "Klebsiella pneumoniae": CriterioQCEspecie(
        especie="Klebsiella pneumoniae",
        tamanho_min_pb=5_000_000,
        tamanho_max_pb=6_200_000,
        max_contigs=500,
        min_n50_pb=20_000,
        max_fracao_n=0.02,
        min_completude=95.0,
        max_contaminacao=5.0,
    ),
    "Staphylococcus aureus": CriterioQCEspecie(
        especie="Staphylococcus aureus",
        tamanho_min_pb=2_600_000,
        tamanho_max_pb=3_200_000,
        max_contigs=300,
        min_n50_pb=25_000,
        max_fracao_n=0.02,
        min_completude=95.0,
        max_contaminacao=5.0,
    ),
}


@dataclass(frozen=True)
class DeterminanteCatalogado:
    """Determinante de resistência catalogado no CARD / AMRFinderPlus."""

    id_atributo: str
    gene: str
    alteracao: str | None
    tipo: str  # "gene_amr" (F1) ou "mutacao_pontual" (F2)
    farmacos_associados: tuple[str, ...]
    fonte: str  # "CARD" ou "AMRFinderPlus"
    mecanismo: str
    sequencia_assinatura: str


# Catálogo curado de determinantes conhecidos (F1 e F2) para E. coli (Fase 1)
CATALOGO_DETERMINANTES_ECOLI: tuple[DeterminanteCatalogado, ...] = (
    DeterminanteCatalogado(
        id_atributo="mut_gyrA_S83L",
        gene="gyrA",
        alteracao="S83L",
        tipo="mutacao_pontual",
        farmacos_associados=("ciprofloxacina",),
        fonte="AMRFinderPlus",
        mecanismo="Alteração do alvo DNA girase subunidade A (QRDR S83L)",
        sequencia_assinatura="CATGGTGACGTATTGTATGCCATGTCGGTCA",
    ),
    DeterminanteCatalogado(
        id_atributo="mut_gyrA_D87N",
        gene="gyrA",
        alteracao="D87N",
        tipo="mutacao_pontual",
        farmacos_associados=("ciprofloxacina",),
        fonte="AMRFinderPlus",
        mecanismo="Alteração do alvo DNA girase subunidade A (QRDR D87N)",
        sequencia_assinatura="GTCATAGTCGGTCGTGCGCTGCCAGACGTGC",
    ),
    DeterminanteCatalogado(
        id_atributo="mut_parC_S80I",
        gene="parC",
        alteracao="S80I",
        tipo="mutacao_pontual",
        farmacos_associados=("ciprofloxacina",),
        fonte="AMRFinderPlus",
        mecanismo="Alteração da topoisomerase IV subunidade A (QRDR S80I)",
        sequencia_assinatura="CCATGGCGATGTCTTATATGCTATGTCGATT",
    ),
    DeterminanteCatalogado(
        id_atributo="gene_qnrS1",
        gene="qnrS1",
        alteracao=None,
        tipo="gene_amr",
        farmacos_associados=("ciprofloxacina",),
        fonte="CARD",
        mecanismo="Proteção da DNA girase mediada por plasmídeo (PMQR)",
        sequencia_assinatura="ATGGAAACCTACAATCATACATATCAACGCA",
    ),
    DeterminanteCatalogado(
        id_atributo="gene_aac6_Ib_cr",
        gene="aac(6')-Ib-cr",
        alteracao=None,
        tipo="gene_amr",
        farmacos_associados=("ciprofloxacina", "gentamicina"),
        fonte="CARD",
        mecanismo="Fluoroquinolona/aminoglicosídeo acetiltransferase bifuncional",
        sequencia_assinatura="ATGACTGAGCATGACCTTGCGATGCTCTATG",
    ),
    DeterminanteCatalogado(
        id_atributo="gene_aac3_IIa",
        gene="aac(3)-IIa",
        alteracao=None,
        tipo="gene_amr",
        farmacos_associados=("gentamicina",),
        fonte="CARD",
        mecanismo="Aminoglicosídeo 3-N-acetiltransferase",
        sequencia_assinatura="ATGCATACGCGGAAGGCAATAACGGAGGCGC",
    ),
    DeterminanteCatalogado(
        id_atributo="gene_aac3_IId",
        gene="aac(3)-IId",
        alteracao=None,
        tipo="gene_amr",
        farmacos_associados=("gentamicina",),
        fonte="AMRFinderPlus",
        mecanismo="Aminoglicosídeo 3-N-acetiltransferase variante IId",
        sequencia_assinatura="ATGCATACGCGGAAGGCAATAACGGAGGCGT",
    ),
    DeterminanteCatalogado(
        id_atributo="gene_ant2_Ia",
        gene="ant(2'')-Ia",
        alteracao=None,
        tipo="gene_amr",
        farmacos_associados=("gentamicina",),
        fonte="CARD",
        mecanismo="Aminoglicosídeo 2''-O-nucleotidiltransferase",
        sequencia_assinatura="ATGGACACAACGCAGGTCGCCTTTGACGCGC",
    ),
    DeterminanteCatalogado(
        id_atributo="gene_armA",
        gene="armA",
        alteracao=None,
        tipo="gene_amr",
        farmacos_associados=("gentamicina",),
        fonte="CARD",
        mecanismo="Metiltransferase de rRNA 16S",
        sequencia_assinatura="ATGGATAAGAATGATGTTGTTAGCATTTTAA",
    ),
    DeterminanteCatalogado(
        id_atributo="gene_blaCTX_M_15",
        gene="blaCTX-M-15",
        alteracao=None,
        tipo="gene_amr",
        farmacos_associados=("cefuroxima", "ceftriaxona"),
        fonte="CARD",
        mecanismo="Beta-lactamase de espectro estendido (ESBL CTX-M-15)",
        sequencia_assinatura="ATGGTTAAAAAATCACTGCGCCAGTTCACGC",
    ),
    DeterminanteCatalogado(
        id_atributo="gene_blaCMY_2",
        gene="blaCMY-2",
        alteracao=None,
        tipo="gene_amr",
        farmacos_associados=("cefuroxima", "ceftriaxona"),
        fonte="AMRFinderPlus",
        mecanismo="Cefalosporinase tipo AmpC plasmidial (CMY-2)",
        sequencia_assinatura="ATGATGAAAAAATCGTTATGCTGCGCTCTGC",
    ),
    DeterminanteCatalogado(
        id_atributo="gene_blaCTX_M_14",
        gene="blaCTX-M-14",
        alteracao=None,
        tipo="gene_amr",
        farmacos_associados=("cefuroxima", "ceftriaxona"),
        fonte="CARD",
        mecanismo="Beta-lactamase de espectro estendido (ESBL CTX-M-14)",
        sequencia_assinatura="ATGGTGACAAAGAGAGTGCAACGGATGATGT",
    ),
    DeterminanteCatalogado(
        id_atributo="mut_ampC_C42T",
        gene="ampC",
        alteracao="c.-42C>T",
        tipo="mutacao_pontual",
        farmacos_associados=("cefuroxima",),
        fonte="AMRFinderPlus",
        mecanismo="Superexpressão de beta-lactamase cromossômica AmpC por mutação promotora",
        sequencia_assinatura="TGACAGTTGTCACGCTGATTGGTGTCGTTAC",
    ),
)

# Catálogo curado para Klebsiella pneumoniae (Fase 2 do escopo de espécies — Seção 3.2)
CATALOGO_DETERMINANTES_KPNEUMONIAE: tuple[DeterminanteCatalogado, ...] = (
    DeterminanteCatalogado(
        id_atributo="mut_gyrA_S83L",
        gene="gyrA",
        alteracao="S83I",
        tipo="mutacao_pontual",
        farmacos_associados=("ciprofloxacina",),
        fonte="AMRFinderPlus",
        mecanismo="Alteração da DNA girase A em K. pneumoniae (QRDR S83I)",
        sequencia_assinatura="CATGGTGACGTATTGTATGCCATGTCGGTCA",
    ),
    DeterminanteCatalogado(
        id_atributo="gene_qnrS1",
        gene="qnrB1",
        alteracao=None,
        tipo="gene_amr",
        farmacos_associados=("ciprofloxacina",),
        fonte="CARD",
        mecanismo="Proteção plasmidial de topoisomerase QnrB1",
        sequencia_assinatura="ATGGAAACCTACAATCATACATATCAACGCA",
    ),
    DeterminanteCatalogado(
        id_atributo="gene_aac3_IIa",
        gene="aac(3)-IIa",
        alteracao=None,
        tipo="gene_amr",
        farmacos_associados=("gentamicina",),
        fonte="CARD",
        mecanismo="Aminoglicosídeo 3-N-acetiltransferase",
        sequencia_assinatura="ATGCATACGCGGAAGGCAATAACGGAGGCGC",
    ),
    DeterminanteCatalogado(
        id_atributo="gene_armA",
        gene="armA",
        alteracao=None,
        tipo="gene_amr",
        farmacos_associados=("gentamicina",),
        fonte="CARD",
        mecanismo="Metiltransferase de rRNA 16S de alto nível",
        sequencia_assinatura="ATGGATAAGAATGATGTTGTTAGCATTTTAA",
    ),
    DeterminanteCatalogado(
        id_atributo="gene_blaCTX_M_15",
        gene="blaCTX-M-15",
        alteracao=None,
        tipo="gene_amr",
        farmacos_associados=("ceftriaxona",),
        fonte="CARD",
        mecanismo="Beta-lactamase de espectro estendido (ESBL CTX-M-15)",
        sequencia_assinatura="ATGGTTAAAAAATCACTGCGCCAGTTCACGC",
    ),
    DeterminanteCatalogado(
        id_atributo="gene_blaCMY_2",
        gene="blaKPC-2",
        alteracao=None,
        tipo="gene_amr",
        farmacos_associados=("ceftriaxona",),
        fonte="CARD",
        mecanismo="Carbapenemase/cefalosporinase classe A KPC-2",
        sequencia_assinatura="ATGATGAAAAAATCGTTATGCTGCGCTCTGC",
    ),
)

# Catálogo curado para Staphylococcus aureus (Fase 3 do escopo de espécies — Seção 3.2: MRSA mecA)
CATALOGO_DETERMINANTES_SAUREUS: tuple[DeterminanteCatalogado, ...] = (
    DeterminanteCatalogado(
        id_atributo="gene_sa_mecA",
        gene="mecA",
        alteracao=None,
        tipo="gene_amr",
        farmacos_associados=("oxacilina",),
        fonte="CARD",
        mecanismo="Proteína ligadora de penicilina PBP2a de baixa afinidade (MRSA)",
        sequencia_assinatura="ATGAAAAAGATAAAAATTGTTCCACTTATTT",
    ),
    DeterminanteCatalogado(
        id_atributo="gene_sa_mecC",
        gene="mecC",
        alteracao=None,
        tipo="gene_amr",
        farmacos_associados=("oxacilina",),
        fonte="AMRFinderPlus",
        mecanismo="Homólogo divergente PBP2a (mecC) em S. aureus",
        sequencia_assinatura="ATGAAAAAAATTTATATTAGTTTAGTTGCGC",
    ),
    DeterminanteCatalogado(
        id_atributo="mut_sa_grlA_S80F",
        gene="grlA",
        alteracao="S80F",
        tipo="mutacao_pontual",
        farmacos_associados=("ciprofloxacina",),
        fonte="AMRFinderPlus",
        mecanismo="Alteração da topoisomerase IV GrlA (QRDR S80F)",
        sequencia_assinatura="CCATGGCGATGTCTTATATGCTATGTCGATT",
    ),
    DeterminanteCatalogado(
        id_atributo="mut_sa_gyrA_S84L",
        gene="gyrA",
        alteracao="S84L",
        tipo="mutacao_pontual",
        farmacos_associados=("ciprofloxacina",),
        fonte="AMRFinderPlus",
        mecanismo="Alteração da DNA girase GyrA em S. aureus (S84L)",
        sequencia_assinatura="CATGGTGACGTATTGTATGCCATGTCGGTCA",
    ),
)

CATALOGO_DETERMINANTES_POR_ESPECIE: dict[str, tuple[DeterminanteCatalogado, ...]] = {
    "Escherichia coli": CATALOGO_DETERMINANTES_ECOLI,
    "Klebsiella pneumoniae": CATALOGO_DETERMINANTES_KPNEUMONIAE,
    "Staphylococcus aureus": CATALOGO_DETERMINANTES_SAUREUS,
}
