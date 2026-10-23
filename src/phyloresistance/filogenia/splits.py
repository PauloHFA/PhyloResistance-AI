"""Estratégias de divisão de dados: Grupo Filogenético, LOLO, Temporal, Geográfica e Aleatória (Seção 7.1, REG-01, REG-09)."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from sklearn.model_selection import StratifiedGroupKFold, StratifiedKFold

from phyloresistance.config import GLOBAL_RANDOM_SEED


@dataclass(frozen=True)
class DivisaoTreinoTeste:
    """Representa uma dobra (fold) ou partição de treino e teste."""

    estrategia: str
    indice_dobra: int
    indices_treino: np.ndarray
    indices_teste: np.ndarray
    grupos_treino: set[str]
    grupos_teste: set[str]
    aviso_vazamento: bool = False
    descricao: str = ""


def verificar_ausencia_vazamento_grupo(
    indices_treino: np.ndarray,
    indices_teste: np.ndarray,
    grupos: list[str] | np.ndarray,
) -> bool:
    """
    Verifica a Regra REG-01: nenhum isolado nem nenhum grupo filogenético pode estar
    simultaneamente no conjunto de treino e no conjunto de teste.
    """
    set_idx_tr = set(map(int, indices_treino))
    set_idx_te = set(map(int, indices_teste))
    if not set_idx_tr.isdisjoint(set_idx_te):
        return False

    arr_grupos = np.asarray(grupos, dtype=str)
    grupos_tr = set(arr_grupos[indices_treino])
    grupos_te = set(arr_grupos[indices_teste])
    return grupos_tr.isdisjoint(grupos_te)


def gerar_splits_grupo_filogenetico(
    y: np.ndarray,
    grupos_filogeneticos: list[str] | np.ndarray,
    n_splits: int = 5,
    seed: int = GLOBAL_RANDOM_SEED,
) -> list[DivisaoTreinoTeste]:
    """
    Estratégia principal de validação (Seção 7.1 e REG-01):
    StratifiedGroupKFold garantindo que todos os membros de um mesmo cluster filogenético/ST
    permaneçam do mesmo lado da divisão.
    """
    y_arr = np.asarray(y, dtype=np.int32)
    grupos_arr = np.asarray(grupos_filogeneticos, dtype=str)
    X_dummy = np.zeros((len(y_arr), 1))

    grupos_unicos = np.unique(grupos_arr)
    n_splits_efetivo = min(n_splits, len(grupos_unicos))
    if n_splits_efetivo < 2:
        raise ValueError("São necessários pelo menos 2 grupos filogenéticos distintos para StratifiedGroupKFold.")

    sgkf = StratifiedGroupKFold(n_splits=n_splits_efetivo, shuffle=True, random_state=seed)
    divisoes: list[DivisaoTreinoTeste] = []

    for idx_dobra, (idx_tr, idx_te) in enumerate(sgkf.split(X_dummy, y_arr, groups=grupos_arr)):
        if not verificar_ausencia_vazamento_grupo(idx_tr, idx_te, grupos_arr):
            raise RuntimeError(f"Violação REG-01 detectada na dobra {idx_dobra}!")

        divisoes.append(
            DivisaoTreinoTeste(
                estrategia="grupo_filogenetico",
                indice_dobra=idx_dobra,
                indices_treino=idx_tr,
                indices_teste=idx_te,
                grupos_treino=set(grupos_arr[idx_tr]),
                grupos_teste=set(grupos_arr[idx_te]),
                aviso_vazamento=False,
                descricao=f"StratifiedGroupKFold (dobra {idx_dobra + 1}/{n_splits_efetivo})",
            )
        )
    return divisoes


def gerar_splits_leave_one_lineage_out(
    y: np.ndarray,
    linhagens_st: list[str] | np.ndarray,
    max_linhagens_avaliadas: int = 5,
) -> list[DivisaoTreinoTeste]:
    """
    Estratégia Leave-One-Lineage-Out (LOLO - Seção 7.1):
    Deixa uma linhagem (ST) inteira de fora em cada partição de teste para avaliar generalização.
    """
    y_arr = np.asarray(y, dtype=np.int32)
    st_arr = np.asarray(linhagens_st, dtype=str)
    valores_st, contagens = np.unique(st_arr, return_counts=True)
    ordem = np.argsort(contagens)[::-1]

    divisoes: list[DivisaoTreinoTeste] = []
    idx_dobra = 0
    for idx_st in ordem:
        st_alvo = str(valores_st[idx_st])
        idx_te = np.where(st_arr == st_alvo)[0]
        idx_tr = np.where(st_arr != st_alvo)[0]
        # Garante que existam ambas as classes no teste e no treino para cálculo de AUROC/AUPRC
        if len(np.unique(y_arr[idx_te])) < 2 or len(np.unique(y_arr[idx_tr])) < 2:
            continue

        divisoes.append(
            DivisaoTreinoTeste(
                estrategia="leave_one_lineage_out",
                indice_dobra=idx_dobra,
                indices_treino=idx_tr,
                indices_teste=idx_te,
                grupos_treino=set(st_arr[idx_tr]),
                grupos_teste={st_alvo},
                aviso_vazamento=False,
                descricao=f"Leave-One-Lineage-Out (teste={st_alvo})",
            )
        )
        idx_dobra += 1
        if idx_dobra >= max_linhagens_avaliadas:
            break

    return divisoes


def gerar_split_temporal(
    anos: list[int] | np.ndarray,
    ano_corte: int = 2022,
) -> DivisaoTreinoTeste:
    """
    Estratégia Temporal (Seção 7.1 e REG-09):
    Treina em isolados coletados até `ano_corte` (inclusive) e testa nos isolados mais recentes (> ano_corte).
    """
    anos_arr = np.asarray(anos, dtype=np.int32)
    idx_tr = np.where(anos_arr <= ano_corte)[0]
    idx_te = np.where(anos_arr > ano_corte)[0]
    if len(idx_tr) == 0 or len(idx_te) == 0:
        mediana = int(np.median(anos_arr))
        idx_tr = np.where(anos_arr <= mediana)[0]
        idx_te = np.where(anos_arr > mediana)[0]
        ano_corte = mediana

    return DivisaoTreinoTeste(
        estrategia="temporal",
        indice_dobra=0,
        indices_treino=idx_tr,
        indices_teste=idx_te,
        grupos_treino={str(a) for a in np.unique(anos_arr[idx_tr])},
        grupos_teste={str(a) for a in np.unique(anos_arr[idx_te])},
        aviso_vazamento=False,
        descricao=f"Treino <= {ano_corte}, Teste > {ano_corte}",
    )


def gerar_split_geografico(
    paises: list[str] | np.ndarray,
    paises_teste: tuple[str, ...] = ("IN", "JP", "FR"),
) -> DivisaoTreinoTeste:
    """
    Estratégia Geográfica (Seção 7.1 e REG-09):
    Treina em um conjunto de países e testa em países retidos para medir robustez à epidemiologia local.
    """
    paises_arr = np.asarray(paises, dtype=str)
    mascara_te = np.isin(paises_arr, paises_teste)
    idx_te = np.where(mascara_te)[0]
    idx_tr = np.where(~mascara_te)[0]

    return DivisaoTreinoTeste(
        estrategia="geografica",
        indice_dobra=0,
        indices_treino=idx_tr,
        indices_teste=idx_te,
        grupos_treino=set(paises_arr[idx_tr]),
        grupos_teste=set(paises_arr[idx_te]),
        aviso_vazamento=False,
        descricao=f"Teste em países retidos: {', '.join(paises_teste)}",
    )


def gerar_splits_aleatorios_vazamento(
    y: np.ndarray,
    grupos_filogeneticos: list[str] | np.ndarray,
    n_splits: int = 5,
    seed: int = GLOBAL_RANDOM_SEED,
) -> list[DivisaoTreinoTeste]:
    """
    Split Aleatório Comum (Seção 7.1):
    Usado APENAS para demonstrar o vazamento filogenético (clones idênticos/mesmo ST em treino e teste)
    e reportar como limite superior enganoso.
    """
    y_arr = np.asarray(y, dtype=np.int32)
    grupos_arr = np.asarray(grupos_filogeneticos, dtype=str)
    X_dummy = np.zeros((len(y_arr), 1))

    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    divisoes: list[DivisaoTreinoTeste] = []
    for idx_dobra, (idx_tr, idx_te) in enumerate(skf.split(X_dummy, y_arr)):
        divisoes.append(
            DivisaoTreinoTeste(
                estrategia="aleatoria_vazamento",
                indice_dobra=idx_dobra,
                indices_treino=idx_tr,
                indices_teste=idx_te,
                grupos_treino=set(grupos_arr[idx_tr]),
                grupos_teste=set(grupos_arr[idx_te]),
                aviso_vazamento=True,
                descricao="LIMITE SUPERIOR ENGANOSO: Split aleatório com vazamento clonal entre treino e teste",
            )
        )
    return divisoes

