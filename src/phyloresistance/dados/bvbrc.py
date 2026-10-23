"""Ingestão de dados do BV-BRC, cálculo de hashes SHA-256 e gerador de coorte piloto realista."""

from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np
import pandas as pd

from phyloresistance.config import (
    CATALOGO_DETERMINANTES_ECOLI,
    CATALOGO_DETERMINANTES_POR_ESPECIE,
    CRITERIOS_QC_POR_ESPECIE,
    ESPECIE_PILOTO,
    FARMACOS_POR_ESPECIE,
    GLOBAL_RANDOM_SEED,
    VERSOES_BANCOS,
)
from phyloresistance.dados.contratos import (
    ManifestoVersaoDados,
    RegistroBrutoBVBRC,
)

# Catálogo de motivos nucleotídicos canônicos (31-mers) reais/sintéticos associados a cada determinante
# Usados na síntese e anotação de sequências FASTA para testes ponta a ponta
ASSINATURAS_NUCLEOTIDICAS_31MER: dict[str, str] = {
    # Determinantes conhecidos E. coli / K. pneumoniae (CARD / AMRFinderPlus)
    "mut_gyrA_S83L": "CATGGTGACGTATTGTATGCCATGTCogtCA".upper().replace("O", "G"),
    "mut_gyrA_D87N": "GTCATAGTCGGTCGTGCGCTGCCAGACGTGC",
    "mut_parC_S80I": "Ccatggcgatgtcttatatgctatgtcgatt".upper(),
    "gene_qnrS1": "ATGGAAACCTACAATCATACATATCAACGCA",
    "gene_aac6_Ib_cr": "ATGACTGAGCATGACCTTGCGATGCTCTATG",
    "gene_aac3_IIa": "ATGCATACGCGGAAGGCAATAACGGAGGCGC",
    "gene_aac3_IId": "ATGCATACGCGGAAGGCAATAACGGAGGCGT",
    "gene_ant2_Ia": "ATGGACACAACGCAGGTCGCCTTTGACGCGC",
    "gene_armA": "ATGGATAAGAATGATGTTGTTAKHATTTTAA".replace("K", "G").replace("H", "C"),
    "gene_blaCTX_M_15": "ATGGTTAAAAAATCACTGCGCCAGTTCACGC",
    "gene_blaCMY_2": "ATGATGAAAAAATCGTTATGCTGCGCTCTGC",
    "gene_blaCTX_M_14": "ATGGTGACAAAGAGAGTGCAACGGATGATGT",
    "mut_ampC_C42T": "TGACAGTTGTCACGCTGATTGGTGTCGTTAC",
    # Determinantes conhecidos S. aureus (MRSA mecA/mecC e grlA/gyrA)
    "gene_sa_mecA": "ATGAAAAAGATAAAAATTGTTCCACTTATTT",
    "gene_sa_mecC": "ATGAAAAAAATTTATATTAGTTTAGTTGCGC",
    "mut_sa_grlA_S80F": "ACATGGCGATGTCTTATATGCTATGTCGATC",
    "mut_sa_gyrA_S84L": "GATGGTGACGTATTGTATGCCATGTCGGTCG",
    # Loci candidatos não catalogados no baseline B1 (capturados por F3 k-mers longos / unitigs)
    "unitig_cand_marR_G103S": "GCGCAAATTCCGGTGCTGGATGAACTGGCGC",
    "unitig_cand_acrR_frameshift": "CGTCTGGCGCTGAAAGAGTTCCGTGACGCGT",
    "unitig_cand_cpxA_Y144C": "ACCCTGGCGCTGGTGITGCTGGCGCTGGCGA".replace("I", "T"),
    "unitig_cand_ompF_prom_del": "TAAACGTTZGTACGTACGTAGCTAGCTAGCA".replace("Z", "C"),
    "unitig_cand_ftsI_insYRIN": "TACCGJATCAACTACCGCATCAACTACCGCA".replace("J", "C"),
    "unitig_cand_sa_pbp4_prom": "GCTAGCTAGCATCGATCGATCGTAGCTAGCG",
    "unitig_cand_sa_norA_prom": "ATCGATCGATCGTAGCTAGCTAGCGATCGAT",
}

MAPEAMENTO_LOCUS_CANDIDATOS: dict[str, dict[str, str]] = {
    "unitig_cand_marR_G103S": {
        "gene": "marR",
        "alteracao": "G103S (unitig)",
        "tipo": "kmer_unitig",
        "fonte": "Descoberta ML (F3)",
        "status": "candidata",
        "mecanismo": "Regulador reprimido de marAB / bomba de efluxo AcrAB-TolC",
    },
    "unitig_cand_acrR_frameshift": {
        "gene": "acrR",
        "alteracao": "fs*45 (unitig)",
        "tipo": "kmer_unitig",
        "fonte": "Descoberta ML (F3)",
        "status": "candidata",
        "mecanismo": "Desrepressão de bomba de efluxo multidroga AcrAB",
    },
    "unitig_cand_cpxA_Y144C": {
        "gene": "cpxA",
        "alteracao": "Y144C (unitig)",
        "tipo": "kmer_unitig",
        "fonte": "Descoberta ML (F3)",
        "status": "candidata",
        "mecanismo": "Sensor quinase de envelope associado a tolerância a aminoglicosídeos",
    },
    "unitig_cand_ompF_prom_del": {
        "gene": "ompF",
        "alteracao": "c.-112del (unitig)",
        "tipo": "kmer_unitig",
        "fonte": "Descoberta ML (F3)",
        "status": "candidata",
        "mecanismo": "Redução de permeabilidade de porina de membrana externa",
    },
    "unitig_cand_ftsI_insYRIN": {
        "gene": "ftsI (PBP3)",
        "alteracao": "ins333YRIN (unitig)",
        "tipo": "kmer_unitig",
        "fonte": "Descoberta ML (F3)",
        "status": "candidata",
        "mecanismo": "Redução de afinidade de PBP3 a cefalosporinas",
    },
    "unitig_cand_sa_pbp4_prom": {
        "gene": "pbp4",
        "alteracao": "c.-90T>A (unitig)",
        "tipo": "kmer_unitig",
        "fonte": "Descoberta ML (F3)",
        "status": "candidata",
        "mecanismo": "Superexpressão de PBP4 conferindo resistência borderline a oxacilina (BORSA)",
    },
    "unitig_cand_sa_norA_prom": {
        "gene": "norA",
        "alteracao": "flqB_up (unitig)",
        "tipo": "kmer_unitig",
        "fonte": "Descoberta ML (F3)",
        "status": "candidata",
        "mecanismo": "Superexpressão da bomba de efluxo de fluoroquinolonas NorA em S. aureus",
    },
}


def calcular_sha256_arquivo(caminho: str | Path) -> str:
    """Calcula o hash SHA-256 de um arquivo em disco para versionamento de dados (Fase 0 / REG-13)."""
    h = hashlib.sha256()
    with open(caminho, "rb") as f:
        for bloco in iter(lambda: f.read(65536), b""):
            h.update(bloco)
    return h.hexdigest()


def gerar_manifesto_dados(caminhos: list[str | Path]) -> ManifestoVersaoDados:
    """Gera manifesto versionado com hashes SHA-256 e versões travadas dos bancos."""
    hashes: dict[str, str] = {}
    for c in caminhos:
        p = Path(c)
        if p.exists() and p.is_file():
            hashes[p.name] = calcular_sha256_arquivo(p)
    return ManifestoVersaoDados(
        versoes_bancos=dict(VERSOES_BANCOS),
        arquivos_hashes_sha256=hashes,
    )


def carregar_tabela_bvbrc(caminho: str | Path) -> list[RegistroBrutoBVBRC]:
    """Carrega registros brutos de fenótipo a partir de arquivo TSV, CSV ou Parquet do BV-BRC."""
    p = Path(caminho)
    if p.suffix.lower() == ".parquet":
        df = pd.read_parquet(p)
    elif p.suffix.lower() in {".tsv", ".txt"}:
        df = pd.read_csv(p, sep="\t")
    else:
        df = pd.read_csv(p)

    registros: list[RegistroBrutoBVBRC] = []
    for row in df.to_dict(orient="records"):
        registros.append(RegistroBrutoBVBRC(**row))
    return registros


def gerar_coorte_piloto_ecoli(
    n_isolados_base: int = 360,
    seed: int = GLOBAL_RANDOM_SEED,
) -> tuple[list[RegistroBrutoBVBRC], pd.DataFrame, pd.DataFrame]:
    """
    Gera uma coorte piloto realista e determinística de Escherichia coli com estrutura clonal (STs),
    distribuição geográfica/temporal, determinantes conhecidos (F1+F2), variantes não catalogadas (F3)
    e casos de borda para testar o funil de exclusão (REG-02, REG-06, REG-07).

    Retorna:
      - registros_brutos: Lista de RegistroBrutoBVBRC (incluindo ruído, conflitos e intermediários).
      - df_genomas_meta: Metadados genômicos e métricas de montagem por isolado para QC/desduplicação.
      - df_genotipos_verdade: Matriz de determinantes F1, F2 e loci candidatos F3 por isolado.
    """
    rng = np.random.default_rng(seed)

    linhagens_st = [
        ("ST131", 0.72),   # Clone pandêmico altamente resistente
        ("ST1193", 0.65),  # Clone emergente resistente a fluoroquinolonas
        ("ST410", 0.60),   # Clone internacional resistente
        ("ST648", 0.55),
        ("ST405", 0.50),
        ("ST38", 0.42),
        ("ST69", 0.28),
        ("ST10", 0.22),
        ("ST73", 0.14),    # Linhagem frequentemente sensível
        ("ST95", 0.10),    # Linhagem frequentemente sensível
        ("ST12", 0.12),
        ("ST127", 0.15),
    ]
    paises = ["BR", "US", "DE", "UK", "FR", "IN", "JP"]
    anos = list(range(2016, 2026))
    metodos = ["Broth microdilution", "Agar dilution", "Disk diffusion"]

    registros_brutos: list[RegistroBrutoBVBRC] = []
    linhas_genomas: list[dict[str, object]] = []
    linhas_genotipos: list[dict[str, int | str]] = []

    ids_det_f1_f2 = [d.id_atributo for d in CATALOGO_DETERMINANTES_ECOLI]
    ids_cand_f3 = list(MAPEAMENTO_LOCUS_CANDIDATOS.keys())

    for i in range(n_isolados_base):
        genome_id = f"ECOLI_{1000 + i:04d}"
        idx_st = i % len(linhagens_st)
        st, prop_res_st = linhagens_st[idx_st]
        pais = paises[int(rng.integers(0, len(paises)))]
        ano = int(rng.choice(anos))
        metodo = str(rng.choice(metodos, p=[0.75, 0.15, 0.10]))

        # Estatísticas de montagem realistas para E. coli (~5.0 Mb)
        tamanho_pb = int(rng.normal(5_050_000, 180_000))
        tamanho_pb = max(4_550_000, min(5_650_000, tamanho_pb))
        contigs = int(rng.integers(25, 190))
        n50 = int(rng.integers(65_000, 320_000))
        fracao_n = float(rng.uniform(0.0001, 0.005))
        completude = float(rng.uniform(97.5, 99.9))
        contaminacao = float(rng.uniform(0.1, 2.2))

        genotipo: dict[str, int | str] = {"genome_id": genome_id, "st": st}
        for col in ids_det_f1_f2 + ids_cand_f3:
            genotipo[col] = 0

        # Sorteio de determinantes condicionado à linhagem ST + variação individual
        # 1) Ciprofloxacina: determinantes conhecidos (F1/F2) vs variantes candidatas (F3)
        tem_cip_conhecido = rng.random() < (prop_res_st * 0.72)
        tem_cip_candidato = (not tem_cip_conhecido) and (rng.random() < (0.20 + 0.15 * prop_res_st))
        if tem_cip_conhecido:
            genotipo["mut_gyrA_S83L"] = 1
            if rng.random() < 0.80:
                genotipo["mut_parC_S80I"] = 1
            if rng.random() < 0.65:
                genotipo["mut_gyrA_D87N"] = 1
            if rng.random() < 0.55:
                genotipo["gene_qnrS1"] = 1
        elif tem_cip_candidato:
            # Resistência mediada por efluxo desreprimido + perda de porina (capturada apenas em F3!)
            genotipo["unitig_cand_marR_G103S"] = 1
            genotipo["unitig_cand_acrR_frameshift"] = 1
            if rng.random() < 0.45:
                genotipo["gene_aac6_Ib_cr"] = 1

        # 2) Gentamicina: determinantes conhecidos (F1) vs variante candidata (F3)
        tem_gen_conhecido = rng.random() < (prop_res_st * 0.68)
        tem_gen_candidato = (not tem_gen_conhecido) and (rng.random() < (0.18 + 0.14 * prop_res_st))
        if tem_gen_conhecido:
            if rng.random() < 0.70:
                genotipo["gene_aac3_IIa"] = 1
            else:
                genotipo["gene_aac3_IId"] = 1
            if rng.random() < 0.55:
                genotipo["gene_ant2_Ia"] = 1
            if rng.random() < 0.35:
                genotipo["gene_armA"] = 1
        elif tem_gen_candidato:
            genotipo["unitig_cand_cpxA_Y144C"] = 1
            genotipo["unitig_cand_acrR_frameshift"] = 1

        # 3) Cefuroxima: determinantes conhecidos (F1/F2) vs variante candidata (F3)
        tem_cxm_conhecido = rng.random() < (prop_res_st * 0.70)
        tem_cxm_candidato = (not tem_cxm_conhecido) and (rng.random() < (0.18 + 0.15 * prop_res_st))
        if tem_cxm_conhecido:
            if rng.random() < 0.70:
                genotipo["gene_blaCTX_M_15"] = 1
            else:
                genotipo["gene_blaCMY_2"] = 1
            if rng.random() < 0.55:
                genotipo["gene_blaCTX_M_14"] = 1
            if rng.random() < 0.45:
                genotipo["mut_ampC_C42T"] = 1
        elif tem_cxm_candidato:
            genotipo["unitig_cand_ompF_prom_del"] = 1
            genotipo["unitig_cand_ftsI_insYRIN"] = 1

        # Determinar fenótipo (MIC e R/S) com penetrância biológica alta (~95%) e baixo ruído basal (~3%)
        for farmaco in ("ciprofloxacina", "gentamicina", "cefuroxima"):
            if farmaco == "ciprofloxacina":
                score_genetico = (
                    2.6 * int(genotipo["mut_gyrA_S83L"])
                    + 2.2 * int(genotipo["mut_gyrA_D87N"])
                    + 2.4 * int(genotipo["mut_parC_S80I"])
                    + 2.0 * int(genotipo["gene_qnrS1"])
                    + 1.5 * int(genotipo["gene_aac6_Ib_cr"])
                    + 2.8 * int(genotipo["unitig_cand_marR_G103S"])
                    + 2.4 * int(genotipo["unitig_cand_acrR_frameshift"])
                )
                prob_r = 1.0 / (1.0 + np.exp(-(score_genetico - 3.2)))
                is_r = bool(rng.random() < prob_r)
                mic = float(rng.choice([2.0, 4.0, 8.0, 16.0])) if is_r else float(rng.choice([0.015, 0.03, 0.06, 0.125]))
            elif farmaco == "gentamicina":
                score_genetico = (
                    3.2 * int(genotipo["gene_aac3_IIa"])
                    + 3.0 * int(genotipo["gene_aac3_IId"])
                    + 2.5 * int(genotipo["gene_ant2_Ia"])
                    + 2.8 * int(genotipo["gene_armA"])
                    + 1.2 * int(genotipo["gene_aac6_Ib_cr"])
                    + 2.8 * int(genotipo["unitig_cand_cpxA_Y144C"])
                    + 2.2 * int(genotipo["unitig_cand_acrR_frameshift"])
                )
                prob_r = 1.0 / (1.0 + np.exp(-(score_genetico - 2.6)))
                is_r = bool(rng.random() < prob_r)
                mic = float(rng.choice([16.0, 32.0, 64.0])) if is_r else float(rng.choice([0.25, 0.5, 1.0]))
            else:  # cefuroxima
                score_genetico = (
                    3.2 * int(genotipo["gene_blaCTX_M_15"])
                    + 3.0 * int(genotipo["gene_blaCMY_2"])
                    + 2.6 * int(genotipo["gene_blaCTX_M_14"])
                    + 2.4 * int(genotipo["mut_ampC_C42T"])
                    + 2.8 * int(genotipo["unitig_cand_ompF_prom_del"])
                    + 2.6 * int(genotipo["unitig_cand_ftsI_insYRIN"])
                )
                prob_r = 1.0 / (1.0 + np.exp(-(score_genetico - 2.6)))
                is_r = bool(rng.random() < prob_r)
                mic = float(rng.choice([32.0, 64.0, 128.0])) if is_r else float(rng.choice([1.0, 2.0, 4.0]))

            registros_brutos.append(
                RegistroBrutoBVBRC(
                    genome_id=genome_id,
                    especie=ESPECIE_PILOTO,
                    farmaco=farmaco,
                    fenotipo_bruto="Resistant" if is_r else "Susceptible",
                    mic_mg_l=mic,
                    metodo_laboratorio=metodo,
                    padrao="CLSI",
                    versao_ponto_corte="2024.1",
                    pais=pais,
                    ano=ano,
                    st=st,
                )
            )

        linhas_genomas.append(
            {
                "genome_id": genome_id,
                "especie_declarada": ESPECIE_PILOTO,
                "especie_inferida": ESPECIE_PILOTO,
                "st": st,
                "pais": pais,
                "ano": ano,
                "tamanho_pb": tamanho_pb,
                "contigs": contigs,
                "n50": n50,
                "fracao_n": fracao_n,
                "completude": completude,
                "contaminacao": contaminacao,
                "clone_de": None,
            }
        )
        linhas_genotipos.append(genotipo)

    # --- Adicionar casos de borda deliberados para exercitar o Funil de Exclusão (REG-02, REG-06, REG-07) ---
    # 1) Espécie fora do escopo (Klebsiella pneumoniae)
    registros_brutos.append(
        RegistroBrutoBVBRC(
            genome_id="KLEB_9001",
            especie="Klebsiella pneumoniae",
            farmaco="ciprofloxacina",
            fenotipo_bruto="Resistant",
            mic_mg_l=4.0,
            pais="BR",
            ano=2023,
            st="ST258",
        )
    )
    # 2) Fármaco fora do escopo (meropenem)
    registros_brutos.append(
        RegistroBrutoBVBRC(
            genome_id="ECOLI_1000",
            especie=ESPECIE_PILOTO,
            farmaco="meropenem",
            fenotipo_bruto="Susceptible",
            mic_mg_l=0.12,
            pais="BR",
            ano=2023,
            st="ST131",
        )
    )
    # 3) Fenótipo Intermediário (I) por MIC (0.5 mg/L para ciprofloxacina)
    registros_brutos.append(
        RegistroBrutoBVBRC(
            genome_id="ECOLI_INT_9002",
            especie=ESPECIE_PILOTO,
            farmaco="ciprofloxacina",
            fenotipo_bruto="Intermediate",
            mic_mg_l=0.5,
            pais="US",
            ano=2023,
            st="ST131",
        )
    )
    # 4) Conflito irreconciliável para um mesmo isolado (uma leitura R e uma S)
    registros_brutos.append(
        RegistroBrutoBVBRC(
            genome_id="ECOLI_CONF_9003",
            especie=ESPECIE_PILOTO,
            farmaco="gentamicina",
            fenotipo_bruto="Resistant",
            mic_mg_l=16.0,
            pais="DE",
            ano=2022,
            st="ST69",
        )
    )
    registros_brutos.append(
        RegistroBrutoBVBRC(
            genome_id="ECOLI_CONF_9003",
            especie=ESPECIE_PILOTO,
            farmaco="gentamicina",
            fenotipo_bruto="Susceptible",
            mic_mg_l=0.5,
            pais="DE",
            ano=2022,
            st="ST69",
        )
    )
    # 5) Genoma reprovado no QC (muitos contigs e baixa completude)
    registros_brutos.append(
        RegistroBrutoBVBRC(
            genome_id="ECOLI_LOWQC_9004",
            especie=ESPECIE_PILOTO,
            farmaco="ciprofloxacina",
            fenotipo_bruto="Resistant",
            mic_mg_l=8.0,
            pais="BR",
            ano=2023,
            st="ST131",
        )
    )
    linhas_genomas.append(
        {
            "genome_id": "ECOLI_LOWQC_9004",
            "especie_declarada": ESPECIE_PILOTO,
            "especie_inferida": ESPECIE_PILOTO,
            "st": "ST131",
            "pais": "BR",
            "ano": 2023,
            "tamanho_pb": 3_100_000,  # abaixo do mínimo de 4.4 Mb
            "contigs": 1250,          # acima do máximo de 500
            "n50": 4_500,             # abaixo de 20.000
            "fracao_n": 0.08,
            "completude": 78.0,
            "contaminacao": 12.5,
            "clone_de": None,
        }
    )
    # 6) Clone quase-idêntico duplicado de ECOLI_1000 para testar desduplicação (REG-02)
    for farmaco in ("ciprofloxacina", "gentamicina", "cefuroxima"):
        reg_orig = next(
            r for r in registros_brutos if r.genome_id == "ECOLI_1000" and r.farmaco == farmaco
        )
        registros_brutos.append(
            RegistroBrutoBVBRC(
                genome_id="ECOLI_DUP_1000",
                especie=ESPECIE_PILOTO,
                farmaco=farmaco,
                fenotipo_bruto=reg_orig.fenotipo_bruto,
                mic_mg_l=reg_orig.mic_mg_l,
                metodo_laboratorio=reg_orig.metodo_laboratorio,
                padrao=reg_orig.padrao,
                versao_ponto_corte=reg_orig.versao_ponto_corte,
                pais=reg_orig.pais,
                ano=reg_orig.ano,
                st=reg_orig.st,
            )
        )
    meta_1000 = dict(linhas_genomas[0])
    meta_1000["genome_id"] = "ECOLI_DUP_1000"
    meta_1000["clone_de"] = "ECOLI_1000"
    linhas_genomas.append(meta_1000)

    geno_1000 = dict(linhas_genotipos[0])
    geno_1000["genome_id"] = "ECOLI_DUP_1000"
    linhas_genotipos.append(geno_1000)

    df_genomas_meta = pd.DataFrame(linhas_genomas)
    df_genotipos_verdade = pd.DataFrame(linhas_genotipos)
    return registros_brutos, df_genomas_meta, df_genotipos_verdade


def gerar_coorte_especie_extensao(
    especie: str,
    n_isolados_base: int = 160,
    seed: int = GLOBAL_RANDOM_SEED,
) -> tuple[list[RegistroBrutoBVBRC], pd.DataFrame, pd.DataFrame]:
    """
    Gera coorte realista para as extensões da Fase 6 (Seção 3.2):
      - Fase 2: Klebsiella pneumoniae (ciprofloxacina, gentamicina, ceftriaxona)
      - Fase 3: Staphylococcus aureus (oxacilina MRSA via mecA, ciprofloxacina via grlA/gyrA)
      - Fase 1: Escherichia coli (delega a gerar_coorte_piloto_ecoli)
    """
    if especie == ESPECIE_PILOTO:
        return gerar_coorte_piloto_ecoli(n_isolados_base=n_isolados_base, seed=seed)

    if especie not in FARMACOS_POR_ESPECIE:
        raise ValueError(f"Espécie '{especie}' não configurada em FARMACOS_POR_ESPECIE.")

    rng = np.random.default_rng(seed + abs(hash(especie)) % 1000)
    farmacos = FARMACOS_POR_ESPECIE[especie]
    catalogo = CATALOGO_DETERMINANTES_POR_ESPECIE[especie]
    crit_qc = CRITERIOS_QC_POR_ESPECIE[especie]

    ids_det_f1_f2 = [d.id_atributo for d in catalogo]
    ids_cand_f3 = list(MAPEAMENTO_LOCUS_CANDIDATOS.keys())

    if especie == "Klebsiella pneumoniae":
        linhagens_st = [
            ("ST258", 0.78),
            ("ST11", 0.75),
            ("ST15", 0.62),
            ("ST147", 0.65),
            ("ST307", 0.58),
            ("ST23", 0.18),
            ("ST37", 0.25),
            ("ST17", 0.20),
        ]
        prefixo = "KPNEU"
    else:  # Staphylococcus aureus
        linhagens_st = [
            ("ST8_USA300", 0.76),
            ("ST5", 0.72),
            ("ST22", 0.68),
            ("ST239", 0.80),
            ("ST30", 0.22),
            ("ST45", 0.18),
            ("ST15", 0.14),
            ("ST1", 0.25),
        ]
        prefixo = "SAUR"

    paises = ["BR", "US", "DE", "UK", "FR", "IN", "JP"]
    anos = list(range(2017, 2026))
    registros_brutos: list[RegistroBrutoBVBRC] = []
    linhas_genomas: list[dict[str, object]] = []
    linhas_genotipos: list[dict[str, int | str]] = []

    tam_medio = (crit_qc.tamanho_min_pb + crit_qc.tamanho_max_pb) // 2

    for i in range(n_isolados_base):
        genome_id = f"{prefixo}_{2000 + i:04d}"
        st, prop_res = linhagens_st[i % len(linhagens_st)]
        pais = paises[int(rng.integers(0, len(paises)))]
        ano = int(rng.choice(anos))

        tamanho_pb = int(
            np.clip(
                rng.normal(tam_medio, 80_000),
                crit_qc.tamanho_min_pb + 50_000,
                crit_qc.tamanho_max_pb - 50_000,
            )
        )
        contigs = int(rng.integers(20, min(150, crit_qc.max_contigs - 20)))
        n50 = int(rng.integers(crit_qc.min_n50_pb + 20_000, 280_000))

        genotipo: dict[str, int | str] = {"genome_id": genome_id, "st": st}
        for col in ids_det_f1_f2 + ids_cand_f3:
            genotipo[col] = 0

        for farmaco in farmacos:
            dets_farm = [d.id_atributo for d in catalogo if farmaco in d.farmacos_associados]
            tem_conhecido = rng.random() < (prop_res * 0.78)
            tem_candidato = (not tem_conhecido) and (rng.random() < (0.16 + 0.12 * prop_res))

            cand_ativo_farmaco: str | None = None
            if tem_conhecido and dets_farm:
                # Ativa o determinante principal (ex.: mecA para oxacilina) e opcionalmente o secundário
                genotipo[dets_farm[0]] = 1
                if len(dets_farm) > 1 and rng.random() < 0.45:
                    genotipo[dets_farm[1]] = 1
            elif tem_candidato:
                if especie == "Staphylococcus aureus":
                    cand_ativo_farmaco = (
                        "unitig_cand_sa_pbp4_prom"
                        if farmaco == "oxacilina"
                        else "unitig_cand_sa_norA_prom"
                    )
                else:
                    if farmaco == "ciprofloxacina":
                        cand_ativo_farmaco = "unitig_cand_marR_G103S"
                        genotipo["unitig_cand_acrR_frameshift"] = 1
                    elif farmaco == "gentamicina":
                        cand_ativo_farmaco = "unitig_cand_cpxA_Y144C"
                    else:
                        cand_ativo_farmaco = "unitig_cand_ompF_prom_del"
                if cand_ativo_farmaco:
                    genotipo[cand_ativo_farmaco] = 1

            soma_conhecidos = sum(int(genotipo[d]) for d in dets_farm)
            tem_cand_farm = 1 if (cand_ativo_farmaco and int(genotipo[cand_ativo_farmaco]) == 1) else 0
            score = 6.0 * soma_conhecidos + 4.8 * tem_cand_farm
            prob_r = 1.0 / (1.0 + np.exp(-(score - 3.2)))
            is_r = bool(rng.random() < prob_r)

            if farmaco == "oxacilina":
                mic = float(rng.choice([8.0, 16.0, 32.0, 64.0])) if is_r else float(rng.choice([0.25, 0.5, 1.0]))
            elif farmaco == "ceftriaxona":
                mic = float(rng.choice([8.0, 16.0, 32.0])) if is_r else float(rng.choice([0.125, 0.25, 0.5]))
            elif farmaco == "gentamicina":
                mic = float(rng.choice([16.0, 32.0, 64.0])) if is_r else float(rng.choice([0.25, 0.5, 1.0]))
            else:
                mic = float(rng.choice([4.0, 8.0, 16.0])) if is_r else float(rng.choice([0.06, 0.125, 0.25]))

            registros_brutos.append(
                RegistroBrutoBVBRC(
                    genome_id=genome_id,
                    especie=especie,
                    farmaco=farmaco,
                    fenotipo_bruto="Resistant" if is_r else "Susceptible",
                    mic_mg_l=mic,
                    metodo_laboratorio="Broth microdilution",
                    padrao="CLSI",
                    versao_ponto_corte="2024.1",
                    pais=pais,
                    ano=ano,
                    st=st,
                )
            )

        linhas_genomas.append(
            {
                "genome_id": genome_id,
                "especie_declarada": especie,
                "especie_inferida": especie,
                "st": st,
                "pais": pais,
                "ano": ano,
                "tamanho_pb": tamanho_pb,
                "contigs": contigs,
                "n50": n50,
                "fracao_n": 0.001,
                "completude": 99.1,
                "contaminacao": 0.6,
                "clone_de": None,
            }
        )
        linhas_genotipos.append(genotipo)

    return registros_brutos, pd.DataFrame(linhas_genomas), pd.DataFrame(linhas_genotipos)


def gerar_coorte_externa_ncbi_pathogen(
    especie: str = ESPECIE_PILOTO,
    n_isolados: int = 120,
    seed: int = 1042,
) -> tuple[list[RegistroBrutoBVBRC], pd.DataFrame, pd.DataFrame]:
    """
    Gera coorte de validação externa independente no formato NCBI Pathogen Detection AST
    (Seção 3.1, 7.1 e Critério de Aceite da Fase 6) com semente e distribuição hospitalar independentes.
    """
    regs, df_meta, df_geno = gerar_coorte_especie_extensao(
        especie=especie, n_isolados_base=n_isolados, seed=seed
    )
    # Filtra apenas isolados base limpos (renomeados com prefixo NCBI_PD_)
    ids_base = set(df_geno["genome_id"].iloc[:n_isolados].values)
    mapa_id = {gid: f"NCBI_PD_{gid}" for gid in ids_base}

    regs_ext: list[RegistroBrutoBVBRC] = []
    for r in regs:
        if r.genome_id in mapa_id and r.especie == especie:
            copia = r.model_copy(update={"genome_id": mapa_id[r.genome_id]})
            regs_ext.append(copia)

    df_meta_ext = df_meta[df_meta["genome_id"].isin(ids_base)].iloc[:n_isolados].copy()
    df_meta_ext["genome_id"] = df_meta_ext["genome_id"].map(mapa_id)

    df_geno_ext = df_geno[df_geno["genome_id"].isin(ids_base)].iloc[:n_isolados].copy()
    df_geno_ext["genome_id"] = df_geno_ext["genome_id"].map(mapa_id)

    return regs_ext, df_meta_ext.reset_index(drop=True), df_geno_ext.reset_index(drop=True)
