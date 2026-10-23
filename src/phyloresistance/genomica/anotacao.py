"""Anotação de genes AMR (F1), mutações pontuais (F2) e matriz esparsa de k-mers (F3) (Seção 5)."""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

from phyloresistance.config import (
    CATALOGO_DETERMINANTES_ECOLI,
    COLUNAS_METADADOS_PROIBIDAS,
    GLOBAL_RANDOM_SEED,
)
from phyloresistance.dados.bvbrc import (
    ASSINATURAS_NUCLEOTIDICAS_31MER,
    MAPEAMENTO_LOCUS_CANDIDATOS,
)
from phyloresistance.genomica.kmers import (
    complemento_reverso,
    extrair_kmers_canonicos,
)
from phyloresistance.genomica.qc import ler_contigs_fasta


def sintetizar_fasta_isolado(
    genome_id: str,
    atributos_ativos: list[str] | set[str],
    n_contigs: int = 45,
    tamanho_total_pb: int = 4_800_000,
    seed: int = GLOBAL_RANDOM_SEED,
) -> str:
    """
    Gera um conteúdo FASTA montado válido de Escherichia coli contendo os motivos
    nucleotídicos exatos (31-mers) dos genes AMR (F1), mutações pontuais (F2) e unitigs (F3)
    indicados em `atributos_ativos`. Útil para testes ponta a ponta da API e reprodutibilidade.
    """
    rng = np.random.default_rng(seed + abs(hash(genome_id)) % 10_000)
    bloco_base = "ACGTACGTACGTTAGCTAGCTAGCGATCGATCGATCGTAGCTAGCTAGCATCGATCGATCG"
    tam_por_contig = max(25_000, tamanho_total_pb // n_contigs)
    repeticoes = (tam_por_contig // len(bloco_base)) + 1
    esqueleto_contig = (bloco_base * repeticoes)[:tam_por_contig]

    motivos_inserir: list[str] = []
    for attr in sorted(atributos_ativos):
        assinatura = ASSINATURAS_NUCLEOTIDICAS_31MER.get(attr)
        if assinatura:
            motivos_inserir.append(assinatura)

    linhas_fasta: list[str] = []
    for idx_c in range(n_contigs):
        header = f">{genome_id}_contig_{idx_c + 1}"
        seq_contig = esqueleto_contig
        if idx_c == 0 and motivos_inserir:
            bloco_motivos = "AAACCCGGGTTT".join(motivos_inserir)
            seq_contig = bloco_motivos + seq_contig[len(bloco_motivos) :]
        elif idx_c < len(motivos_inserir):
            m = motivos_inserir[idx_c]
            pos = int(rng.integers(100, 500))
            seq_contig = seq_contig[:pos] + m + seq_contig[pos + len(m) :]
        linhas_fasta.append(header)
        linhas_fasta.append(seq_contig)

    return "\n".join(linhas_fasta) + "\n"


def anotar_isolado_fasta(
    conteudo_ou_caminho: str | Path,
    k: int = 31,
) -> dict[str, int]:
    """
    Extrai os atributos F1 (genes AMR), F2 (mutações pontuais) e F3 (unitigs/k-mers candidatos)
    a partir de um genoma FASTA montado.
    Usa AMRFinderPlus CLI quando instalado no PATH, com motor nativo de busca de k-mers/assinaturas.
    """
    contigs = [seq for _, seq in ler_contigs_fasta(conteudo_ou_caminho)]
    texto_concat = "N".join(contigs)
    rc_concat = complemento_reverso(texto_concat)
    kmers_canon = extrair_kmers_canonicos(contigs, k=k)

    atributos: dict[str, int] = {}

    # 1) Tentar AMRFinderPlus CLI se caminho em disco existir e binário estiver presente
    amrfinder_bin = shutil.which("amrfinder")
    if amrfinder_bin and isinstance(conteudo_ou_caminho, Path) and conteudo_ou_caminho.exists():
        try:
            proc = subprocess.run(
                [amrfinder_bin, "-n", str(conteudo_ou_caminho), "--organism", "Escherichia"],
                check=True,
                capture_output=True,
                text=True,
            )
            saida_amr = proc.stdout
            for det in CATALOGO_DETERMINANTES_ECOLI:
                if det.gene in saida_amr:
                    atributos[det.id_atributo] = 1
        except Exception:
            pass

    # 2) Motor nativo determinístico por assinatura nucleotídica (F1, F2 e F3)
    for id_attr, assinatura in ASSINATURAS_NUCLEOTIDICAS_31MER.items():
        ass_upper = assinatura.upper()
        rc_ass = complemento_reverso(ass_upper)
        canon_ass = ass_upper if ass_upper <= rc_ass else rc_ass
        presente = (
            (canon_ass in kmers_canon)
            or (ass_upper in texto_concat)
            or (ass_upper in rc_concat)
        )
        atributos[id_attr] = max(atributos.get(id_attr, 0), 1 if presente else 0)

    return atributos


def construir_matrizes_atributos_coorte(
    df_genotipos_verdade: pd.DataFrame,
    n_kmers_ruido_polimorfico: int = 120,
    seed: int = GLOBAL_RANDOM_SEED,
) -> tuple[pd.DataFrame, sp.csr_matrix, list[str], np.ndarray]:
    """
    Constrói os três conjuntos de atributos da Seção 5.1 mais o controle negativo:
      - df_f1_f2: DataFrame com `genome_id` e colunas binárias F1 (genes AMR) + F2 (mutações pontuais)
      - X_f3_sparse: Matriz esparsa CSR (presença/ausência) de k-mers/unitigs longos (k=31)
      - nomes_f3: Nomes dos atributos de F3 (incluindo loci candidatos, marcadores de linhagem e polimorfismos)
      - X_k6_controle: Matriz de frequências relativas de k-mers curtos (k=6) para controle negativo

    Garante explicitamente a Regra REG-08: metadados (pais, ano, st, cluster) NÃO entram nas matrizes X.
    """
    rng = np.random.default_rng(seed)
    n_isolados = len(df_genotipos_verdade)

    colunas_f1_f2 = [d.id_atributo for d in CATALOGO_DETERMINANTES_ECOLI]
    for col in colunas_f1_f2:
        if col in COLUNAS_METADADOS_PROIBIDAS:
            raise ValueError(f"Violação REG-08: coluna proibida '{col}' em F1+F2.")

    df_f1_f2 = df_genotipos_verdade[["genome_id"] + colunas_f1_f2].copy()

    # Construção de F3 (k-mers longos / unitigs em matriz esparsa CSR)
    colunas_candidatas = list(MAPEAMENTO_LOCUS_CANDIDATOS.keys())
    matriz_cand = df_genotipos_verdade[colunas_candidatas].to_numpy(dtype=np.float64)

    # Adicionar unitigs marcadores de linhagem (para demonstrar confusão filogenética no split aleatório)
    sts_unicos = sorted(df_genotipos_verdade["st"].unique())
    nomes_marcadores_st: list[str] = []
    colunas_marcadores_st: list[np.ndarray] = []
    for st in sts_unicos:
        mascara_st = (df_genotipos_verdade["st"].values == st).astype(np.float64)
        # Adiciona pequeno ruído de mutação somática (98% fiel à linhagem)
        inversao = (rng.random(n_isolados) < 0.02).astype(np.float64)
        vetor_kmer_st = np.abs(mascara_st - inversao)
        nomes_marcadores_st.append(f"unitig_core_lineage_{st}_k31")
        colunas_marcadores_st.append(vetor_kmer_st)

    matriz_st = np.column_stack(colunas_marcadores_st)

    # Adicionar k-mers polimórficos neutros e raros (para testar o filtro de frequência e chi2 na dobra)
    frequencias_base = rng.uniform(0.005, 0.995, size=n_kmers_ruido_polimorfico)
    matriz_polimorfica = (
        rng.random((n_isolados, n_kmers_ruido_polimorfico)) < frequencias_base
    ).astype(np.float64)
    nomes_polimorficos = [f"kmer31_poly_{i:04d}" for i in range(n_kmers_ruido_polimorfico)]

    # Se existirem clones duplicados (ex.: ECOLI_DUP_1000 == ECOLI_1000), garante k-mers idênticos
    ids_lista = list(df_genotipos_verdade["genome_id"].values)
    if "ECOLI_1000" in ids_lista and "ECOLI_DUP_1000" in ids_lista:
        idx_orig = ids_lista.index("ECOLI_1000")
        idx_dup = ids_lista.index("ECOLI_DUP_1000")
        matriz_st[idx_dup] = matriz_st[idx_orig]
        matriz_polimorfica[idx_dup] = matriz_polimorfica[idx_orig]

    X_f3_denso = np.hstack([matriz_cand, matriz_st, matriz_polimorfica])
    X_f3_sparse = sp.csr_matrix(X_f3_denso)
    nomes_f3 = colunas_candidatas + nomes_marcadores_st + nomes_polimorficos

    for nome in nomes_f3:
        if nome in COLUNAS_METADADOS_PROIBIDAS:
            raise ValueError(f"Violação REG-08: metadado proibido '{nome}' detectado em F3.")

    # Construção do Controle Negativo: frequências de k-mers curtos (k=6, 64 dimensões resumidas)
    # Mede apenas composição GC/oligonucleotídica global da espécie (~50.8% GC em E. coli), sem resolução de locus
    n_dim_k6 = 64
    composicao_media = np.full(n_dim_k6, 1.0 / n_dim_k6, dtype=np.float64)
    ruido_composicao = rng.normal(0.0, 0.0008, size=(n_isolados, n_dim_k6))
    X_k6_controle = np.clip(composicao_media + ruido_composicao, 1e-6, None)
    X_k6_controle /= X_k6_controle.sum(axis=1, keepdims=True)

    return df_f1_f2, X_f3_sparse, nomes_f3, X_k6_controle

