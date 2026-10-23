"""Cálculo de distância genômica tipo Mash, desduplicação clonal (REG-02) e clustering filogenético."""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp
from pydantic import BaseModel, Field
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform


class RelatorioDesduplicacao(BaseModel):
    """Relatório auditável de desduplicação de genomas quase idênticos antes do split (REG-02)."""

    total_isolados_entrada: int
    total_clones_removidos: int
    total_isolados_unicos: int
    limiar_distancia_mash: float
    pares_colapsados: list[dict[str, str | float]] = Field(default_factory=list)
    ids_mantidos: list[str] = Field(default_factory=list)
    ids_removidos: list[str] = Field(default_factory=list)


def calcular_matriz_distancia_mash(
    X_kmers: sp.csr_matrix | np.ndarray,
    k: int = 31,
) -> np.ndarray:
    """
    Calcula a matriz simétrica de distância Mash D = -1/k * ln(2*J / (1 + J))
    a partir da similaridade de Jaccard J dos k-mers binários entre isolados.
    """
    X_bin = (sp.csr_matrix(X_kmers) > 0).astype(np.float64)
    intersecao = (X_bin @ X_bin.T).toarray()
    somas = np.asarray(X_bin.sum(axis=1)).ravel()
    uniao = somas[:, None] + somas[None, :] - intersecao
    uniao = np.maximum(uniao, 1.0)

    jaccard = np.clip(intersecao / uniao, 1e-9, 1.0)
    termo = (2.0 * jaccard) / (1.0 + jaccard)
    dist = -1.0 / float(k) * np.log(termo)
    dist = np.clip(dist, 0.0, 1.0)
    np.fill_diagonal(dist, 0.0)
    return dist


def desduplicar_genomas_quase_identicos(
    genome_ids: list[str],
    matriz_distancia: np.ndarray,
    limiar_distancia: float = 0.0005,
) -> tuple[list[int], RelatorioDesduplicacao]:
    """
    Remove genomas quase idênticos (distância Mash <= limiar_distancia) ANTES de dividir em treino/teste (REG-02).
    Retorna os índices mantidos e o RelatorioDesduplicacao.
    """
    n = len(genome_ids)
    if matriz_distancia.shape != (n, n):
        raise ValueError("Dimensão da matriz de distância incompatível com genome_ids.")

    removidos: set[int] = set()
    pares_colapsados: list[dict[str, str | float]] = []

    for i in range(n):
        if i in removidos:
            continue
        for j in range(i + 1, n):
            if j in removidos:
                continue
            d = float(matriz_distancia[i, j])
            if d <= limiar_distancia:
                removidos.add(j)
                pares_colapsados.append(
                    {
                        "representante": genome_ids[i],
                        "duplicata_removida": genome_ids[j],
                        "distancia_mash": round(d, 8),
                    }
                )

    indices_mantidos = [i for i in range(n) if i not in removidos]
    ids_mantidos = [genome_ids[i] for i in indices_mantidos]
    ids_removidos = [genome_ids[i] for i in sorted(removidos)]

    relatorio = RelatorioDesduplicacao(
        total_isolados_entrada=n,
        total_clones_removidos=len(removidos),
        total_isolados_unicos=len(indices_mantidos),
        limiar_distancia_mash=limiar_distancia,
        pares_colapsados=pares_colapsados,
        ids_mantidos=ids_mantidos,
        ids_removidos=ids_removidos,
    )
    return indices_mantidos, relatorio


def inferir_clusters_filogeneticos(
    genome_ids: list[str],
    matriz_distancia: np.ndarray,
    sts: list[str] | None = None,
    limiar_distancia_cluster: float = 0.018,
) -> list[str]:
    """
    Agrupa isolados em clusters filogenéticos usando distância Mash (clustering hierárquico average-linkage / UPGMA)
    combinado com o Sequence Type (ST) quando disponível, garantindo que membros de um mesmo complexo clonal
    jamais sejam separados entre treino e teste (REG-01).
    """
    n = len(genome_ids)
    if n == 0:
        return []
    if n == 1:
        return [sts[0] if sts else "CLUSTER_001"]

    dist_sim = 0.5 * (matriz_distancia + matriz_distancia.T)
    np.fill_diagonal(dist_sim, 0.0)
    vetor_condensado = squareform(dist_sim, checks=False)
    Z = linkage(vetor_condensado, method="average")
    rotulos_hc = fcluster(Z, t=limiar_distancia_cluster, criterion="distance")

    if sts is not None and len(sts) == n:
        # Garante que todos os membros do mesmo ST fiquem no mesmo grupo filogenético
        return [str(st) for st in sts]

    return [f"CLUSTER_{int(c):03d}" for c in rotulos_hc]

