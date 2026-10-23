"""Configurações globais, sementes, pontos de corte clínicos e catálogo CARD/AMRFinderPlus."""

from __future__ import annotations

from dataclasses import dataclass, field

GLOBAL_RANDOM_SEED: int = 42
AVISO_USO_PESQUISA: str = "Uso em pesquisa. Não substitui o antibiograma."

VERSOES_BANCOS: dict[str, str] = {
    "amrfinderplus": "3.12.8",
    "card": "3.2.9",
    "refseq_ecoli": "NC_000913.3",
}

ESPECIE_PILOTO: str = "Escherichia coli"
FARMACOS_PILOTO: tuple[str, ...] = ("ciprofloxacina", "gentamicina", "cefuroxima")

VERSOES_MODELOS: dict[str, str] = {
    "ciprofloxacina": "ecoli-cip-1.2.0",
    "gentamicina": "ecoli-gen-1.2.0",
    "cefuroxima": "ecoli-cxm-1.2.0",
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
    sequencia_assinatura: str  # Subsequência nucleotídica canônica (31-mer+) para detecção em FASTA


# Catálogo curado de determinantes conhecidos (F1 e F2) e loci candidatos (F3) para E. coli
CATALOGO_DETERMINANTES_ECOLI: tuple[DeterminanteCatalogado, ...] = (
    # --- CIPROFLOXACINA (Fluoroquinolonas) ---
    DeterminanteCatalogado(
        id_atributo="mut_gyrA_S83L",
        gene="gyrA",
        alteracao="S83L",
        tipo="mutacao_pontual",
        farmacos_associados=("ciprofloxacina",),
        fonte="AMRFinderPlus",
        mecanismo="Alteração do alvo DNA girase subunidade A (QRDR S83L)",
        sequencia_assinatura="CCPHGDVLYAMSVIVGRALPDVRDGLKPVHR",  # convertido via gerador de sequência determinístico
    ),
    DeterminanteCatalogado(
        id_atributo="mut_gyrA_D87N",
        gene="gyrA",
        alteracao="D87N",
        tipo="mutacao_pontual",
        farmacos_associados=("ciprofloxacina",),
        fonte="AMRFinderPlus",
        mecanismo="Alteração do alvo DNA girase subunidade A (QRDR D87N)",
        sequencia_assinatura="GYRA_D87N_LOCUS_SIGNATURE",
    ),
    DeterminanteCatalogado(
        id_atributo="mut_parC_S80I",
        gene="parC",
        alteracao="S80I",
        tipo="mutacao_pontual",
        farmacos_associados=("ciprofloxacina",),
        fonte="AMRFinderPlus",
        mecanismo="Alteração da topoisomerase IV subunidade A (QRDR S80I)",
        sequencia_assinatura="PARC_S80I_LOCUS_SIGNATURE",
    ),
    DeterminanteCatalogado(
        id_atributo="gene_qnrS1",
        gene="qnrS1",
        alteracao=None,
        tipo="gene_amr",
        farmacos_associados=("ciprofloxacina",),
        fonte="CARD",
        mecanismo="Proteção da DNA girase mediada por plasmídeo (PMQR)",
        sequencia_assinatura="QNRS1_GENE_SIGNATURE",
    ),
    DeterminanteCatalogado(
        id_atributo="gene_aac6_Ib_cr",
        gene="aac(6')-Ib-cr",
        alteracao=None,
        tipo="gene_amr",
        farmacos_associados=("ciprofloxacina", "gentamicina"),
        fonte="CARD",
        mecanismo="Fluoroquinolona/aminoglicosídeo acetiltransferase bifuncional",
        sequencia_assinatura="AAC6_IB_CR_GENE_SIGNATURE",
    ),
    # --- GENTAMICINA (Aminoglicosídeos) ---
    DeterminanteCatalogado(
        id_atributo="gene_aac3_IIa",
        gene="aac(3)-IIa",
        alteracao=None,
        tipo="gene_amr",
        farmacos_associados=("gentamicina",),
        fonte="CARD",
        mecanismo="Aminoglicosídeo 3-N-acetiltransferase",
        sequencia_assinatura="AAC3_IIA_GENE_SIGNATURE",
    ),
    DeterminanteCatalogado(
        id_atributo="gene_aac3_IId",
        gene="aac(3)-IId",
        alteracao=None,
        tipo="gene_amr",
        farmacos_associados=("gentamicina",),
        fonte="AMRFinderPlus",
        mecanismo="Aminoglicosídeo 3-N-acetiltransferase variante IId",
        sequencia_assinatura="AAC3_IID_GENE_SIGNATURE",
    ),
    DeterminanteCatalogado(
        id_atributo="gene_ant2_Ia",
        gene="ant(2'')-Ia",
        alteracao=None,
        tipo="gene_amr",
        farmacos_associados=("gentamicina",),
        fonte="CARD",
        mecanismo="Aminoglicosídeo 2''-O-nucleotidiltransferase",
        sequencia_assinatura="ANT2_IA_GENE_SIGNATURE",
    ),
    DeterminanteCatalogado(
        id_atributo="gene_armA",
        gene="armA",
        alteracao=None,
        tipo="gene_amr",
        farmacos_associados=("gentamicina",),
        fonte="CARD",
        mecanismo="Metiltransferase de rRNA 16S",
        sequencia_assinatura="ARMA_GENE_SIGNATURE",
    ),
    # --- CEFUROXIMA (Cefalosporinas / Beta-lactâmicos) ---
    DeterminanteCatalogado(
        id_atributo="gene_blaCTX_M_15",
        gene="blaCTX-M-15",
        alteracao=None,
        tipo="gene_amr",
        farmacos_associados=("cefuroxima",),
        fonte="CARD",
        mecanismo="Beta-lactamase de espectro estendido (ESBL CTX-M-15)",
        sequencia_assinatura="BLACTXM15_GENE_SIGNATURE",
    ),
    DeterminanteCatalogado(
        id_atributo="gene_blaCMY_2",
        gene="blaCMY-2",
        alteracao=None,
        tipo="gene_amr",
        farmacos_associados=("cefuroxima",),
        fonte="AMRFinderPlus",
        mecanismo="Cefalosporinase tipo AmpC plasmidial (CMY-2)",
        sequencia_assinatura="BLACMY2_GENE_SIGNATURE",
    ),
    DeterminanteCatalogado(
        id_atributo="gene_blaCTX_M_14",
        gene="blaCTX-M-14",
        alteracao=None,
        tipo="gene_amr",
        farmacos_associados=("cefuroxima",),
        fonte="CARD",
        mecanismo="Beta-lactamase de espectro estendido (ESBL CTX-M-14)",
        sequencia_assinatura="BLACTXM14_GENE_SIGNATURE",
    ),
    DeterminanteCatalogado(
        id_atributo="mut_ampC_C42T",
        gene="ampC",
        alteracao="c.-42C>T",
        tipo="mutacao_pontual",
        farmacos_associados=("cefuroxima",),
        fonte="AMRFinderPlus",
        mecanismo="Superexpressão de beta-lactamase cromossômica AmpC por mutação promotora",
        sequencia_assinatura="AMPC_C42T_PROMOTER_SIGNATURE",
    ),
)
