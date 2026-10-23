"""Extração de k-mers canônicos (k=21..31), agrupamento em unitigs e seleção anti-vazamento na dobra (Seção 5.2, REG-03)."""

from __future__ import annotations

import hashlib
import shutil
import subprocess
from itertools import product
from pathlib import Path

import numpy as np
import scipy.sparse as sp
from sklearn.feature_selection import chi2

TABELA_COMPLEMENTO = str.maketrans({"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"})


def complemento_reverso(seq: str) -> str:
    """Retorna o complemento reverso de uma sequência de DNA."""
    return seq.upper().translate(TABELA_COMPLEMENTO)[::-1]


def kmer_canonico(kmer: str) -> str:
    """Retorna a forma canônica lexicográfica min(kmer, complemento_reverso(kmer))."""
    rc = complemento_reverso(kmer)
    return kmer if kmer <= rc else rc


def extrair_kmers_canonicos(
    sequencia_ou_contigs: str | list[str],
    k: int = 31,
) -> set[str]:
    """
    Extrai o conjunto binário (presença/ausência) de k-mers canônicos longos (21 <= k <= 31),
    ignorando janelas que contenham bases ambíguas 'N' (Seção 5.2).
    """
    if k < 1:
        raise ValueError("Tamanho k deve ser positivo.")

    contigs = [sequencia_ou_contigs] if isinstance(sequencia_ou_contigs, str) else sequencia_ou_contigs
    kmers_presentes: set[str] = set()

    for contig in contigs:
        seq = contig.upper()
        n = len(seq)
        if n < k:
            continue
        # Divide em blocos livres de 'N' para acelerar a varredura
        for subseq in seq.split("N"):
            tam_sub = len(subseq)
            if tam_sub < k:
                continue
            rc_sub = complemento_reverso(subseq)
            for i in range(tam_sub - k + 1):
                foward_k = subseq[i : i + k]
                rev_k = rc_sub[tam_sub - i - k : tam_sub - i]
                kmers_presentes.add(foward_k if foward_k <= rev_k else rev_k)

    return kmers_presentes


def extrair_kmers_via_kmc_ou_nativo(
    caminho_fasta: str | Path,
    k: int = 31,
    diretorio_tmp: str | Path | None = None,
) -> set[str]:
    """
    Wrapper híbrido: utiliza KMC3 ('kmc' + 'kmc_tools') quando disponível no sistema,
    com fallback automático para o extrator nativo em Python.
    """
    caminho = Path(caminho_fasta)
    kmc_bin = shutil.which("kmc")
    kmc_tools_bin = shutil.which("kmc_tools")

    if kmc_bin and kmc_tools_bin and diretorio_tmp is not None:
        tmp_dir = Path(diretorio_tmp)
        tmp_dir.mkdir(parents=True, exist_ok=True)
        prefixo_saida = tmp_dir / f"kmc_{caminho.stem}"
        arquivo_dump = tmp_dir / f"kmc_{caminho.stem}.txt"
        try:
            subprocess.run(
                [kmc_bin, f"-k{k}", "-ci1", "-fm", "-t4", str(caminho), str(prefixo_saida), str(tmp_dir)],
                check=True,
                capture_output=True,
                text=True,
            )
            subprocess.run(
                [kmc_tools_bin, "transform", str(prefixo_saida), "dump", str(arquivo_dump)],
                check=True,
                capture_output=True,
                text=True,
            )
            kmers: set[str] = set()
            for linha in arquivo_dump.read_text(encoding="utf-8").splitlines():
                partes = linha.strip().split()
                if partes:
                    kmers.add(partes[0].upper())
            return kmers
        except Exception:
            pass

    from phyloresistance.genomica.qc import ler_contigs_fasta

    contigs = [seq for _, seq in ler_contigs_fasta(caminho)]
    return extrair_kmers_canonicos(contigs, k=k)


def extrair_frequencias_kmers_curtos(
    sequencia_ou_contigs: str | list[str],
    k: int = 6,
) -> np.ndarray:
    """
    Extrai vetor de frequências relativas de k-mers curtos (k=6 a 8).
    Usado ESTRITAMENTE como controle negativo (ablação da Seção 2 e Seção 6.1)
    para demonstrar que composição global não captura mutações pontuais.
    """
    contigs = [sequencia_ou_contigs] if isinstance(sequencia_ou_contigs, str) else sequencia_ou_contigs
    todos_kmers = ["".join(p) for p in product("ACGT", repeat=k)]
    indice_kmer = {km: idx for idx, km in enumerate(todos_kmers)}
    contagens = np.zeros(len(todos_kmers), dtype=np.float64)

    for contig in contigs:
        seq = contig.upper()
        for subseq in seq.split("N"):
            for i in range(len(subseq) - k + 1):
                km = subseq[i : i + k]
                idx = indice_kmer.get(km)
                if idx is not None:
                    contagens[idx] += 1.0

    total = contagens.sum()
    if total > 0:
        contagens /= total
    return contagens


def agrupar_kmers_equivalentes_unitigs(
    X_kmers: sp.csr_matrix | np.ndarray,
    nomes_kmers: list[str],
) -> tuple[sp.csr_matrix, list[str], dict[str, list[str]]]:
    """
    Agrupa k-mers equivalentes (que têm exatamente o mesmo padrão de presença/ausência
    nos isolados de referência) em unidades 'unitigs' / padrões (Seção 5.2).
    Reduz dimensão e evita divisão arbitrária de importância SHAP entre k-mers co-ocorrentes.
    """
    matriz_csc = sp.csc_matrix(X_kmers)
    assinaturas: dict[bytes, int] = {}
    colunas_unicas: list[int] = []
    nomes_unitigs: list[str] = []
    mapa_unitig_para_kmers: dict[str, list[str]] = {}

    for j in range(matriz_csc.shape[1]):
        inicio, fim = matriz_csc.indptr[j], matriz_csc.indptr[j + 1]
        indices_linhas = matriz_csc.indices[inicio:fim]
        chave = indices_linhas.tobytes()
        nome_k = nomes_kmers[j]

        if chave not in assinaturas:
            idx_unitig = len(colunas_unicas)
            assinaturas[chave] = idx_unitig
            colunas_unicas.append(j)
            # Se o nome já for semântico (ex.: unitig_cand_*), preserva-o
            nome_u = nome_k if nome_k.startswith("unitig_") else f"unitig_{idx_unitig:04d}_{nome_k[:12]}"
            nomes_unitigs.append(nome_u)
            mapa_unitig_para_kmers[nome_u] = [nome_k]
        else:
            idx_existente = assinaturas[chave]
            nome_u_existente = nomes_unitigs[idx_existente]
            mapa_unitig_para_kmers[nome_u_existente].append(nome_k)

    X_reduzido = matriz_csc[:, colunas_unicas].tocsr()
    return X_reduzido, nomes_unitigs, mapa_unitig_para_kmers


class FiltroSelecaoKmersDobra:
    """
    Transformador compatível com scikit-learn que aplica:
    1) Filtro de frequência mínima e máxima de presença/ausência
    2) Redução supervisionada por teste qui-quadrado (chi2)
    ESTRITAMENTE ajustado nos dados de treino da dobra (REG-03).
    """

    def __init__(
        self,
        min_freq: float = 0.02,
        max_freq: float = 0.98,
        max_kmers_selecionados: int = 150,
    ) -> None:
        self.min_freq = min_freq
        self.max_freq = max_freq
        self.max_kmers_selecionados = max_kmers_selecionados
        self.mascara_frequencia_: np.ndarray | None = None
        self.indices_selecionados_: np.ndarray | None = None
        self.scores_chi2_: np.ndarray | None = None
        self.fingerprint_treino_: str | None = None
        self.n_amostras_treino_: int = 0

    def fit(
        self,
        X_treino: sp.csr_matrix | np.ndarray,
        y_treino: np.ndarray,
        ids_treino: list[str] | None = None,
    ) -> FiltroSelecaoKmersDobra:
        """Ajusta o filtro de frequência e seleção supervisionada APENAS nos dados de treino."""
        X_csr = sp.csr_matrix(X_treino)
        y_arr = np.asarray(y_treino, dtype=np.int32)
        n_amostras = X_csr.shape[0]
        if n_amostras == 0:
            raise ValueError("X_treino não pode ser vazio.")

        self.n_amostras_treino_ = n_amostras
        if ids_treino is not None:
            concat = "|".join(sorted(ids_treino))
            self.fingerprint_treino_ = hashlib.sha256(concat.encode("utf-8")).hexdigest()
        else:
            self.fingerprint_treino_ = hashlib.sha256(y_arr.tobytes()).hexdigest()

        # Presença/ausência binarizada
        X_bin = (X_csr > 0).astype(np.float64)
        freq = np.asarray(X_bin.mean(axis=0)).ravel()

        self.mascara_frequencia_ = (freq >= self.min_freq) & (freq <= self.max_freq)
        indices_freq = np.where(self.mascara_frequencia_)[0]

        if len(indices_freq) == 0:
            # Fallback seguro caso todos os atributos sejam raros em um mini-teste
            indices_freq = np.arange(X_csr.shape[1])
            self.mascara_frequencia_ = np.ones(X_csr.shape[1], dtype=bool)

        X_filtrado = X_bin[:, indices_freq]
        if X_filtrado.shape[1] <= self.max_kmers_selecionados:
            self.indices_selecionados_ = indices_freq
            self.scores_chi2_ = np.ones(len(indices_freq), dtype=np.float64)
        else:
            estat_chi2, _ = chi2(X_filtrado, y_arr)
            estat_chi2 = np.nan_to_num(estat_chi2, nan=0.0)
            ordem = np.argsort(estat_chi2)[::-1][: self.max_kmers_selecionados]
            ordem_ordenada = np.sort(ordem)
            self.indices_selecionados_ = indices_freq[ordem_ordenada]
            self.scores_chi2_ = estat_chi2[ordem_ordenada]

        return self

    def transform(self, X: sp.csr_matrix | np.ndarray) -> sp.csr_matrix:
        """Aplica a máscara de atributos selecionados na dobra de treino sem reestimar estatísticas."""
        if self.indices_selecionados_ is None:
            raise RuntimeError("FiltroSelecaoKmersDobra precisa ser ajustado com fit() antes de transform().")
        X_csr = sp.csr_matrix(X)
        return (X_csr[:, self.indices_selecionados_] > 0).astype(np.float64).tocsr()

    def fit_transform(
        self,
        X_treino: sp.csr_matrix | np.ndarray,
        y_treino: np.ndarray,
        ids_treino: list[str] | None = None,
    ) -> sp.csr_matrix:
        return self.fit(X_treino, y_treino, ids_treino=ids_treino).transform(X_treino)

    def filtrar_nomes(self, nomes_atributos: list[str]) -> list[str]:
        if self.indices_selecionados_ is None:
            raise RuntimeError("FiltroSelecaoKmersDobra não ajustado.")
        return [nomes_atributos[int(i)] for i in self.indices_selecionados_]

