"""Métricas de avaliação (AUROC, AUPRC, VME, ME, Brier, Cobertura Conformal) e Bootstrap em Blocos (Seção 7.2, REG-04, REG-05)."""

from __future__ import annotations

import numpy as np
from pydantic import BaseModel
from sklearn.metrics import (
    average_precision_score,
    balanced_accuracy_score,
    brier_score_loss,
    roc_auc_score,
)

from phyloresistance.config import GLOBAL_RANDOM_SEED
from phyloresistance.modelos.calibracao import PredicaoComIncerteza


class IntervaloConfiancaBloco(BaseModel):
    """Estimativa pontual e intervalo de confiança de 95% via bootstrap em blocos (clusters)."""

    estimativa: float
    ic95_inferior: float
    ic95_superior: float


class ComparacaoPareadaBaseline(BaseModel):
    """Resultado da comparação pareada por cluster entre um modelo de ML e o Baseline B1 (REG-04)."""

    nome_modelo: str
    delta_auprc_medio: float
    delta_auprc_ic95_inf: float
    delta_auprc_ic95_sup: float
    delta_vme_medio: float
    supera_b1_com_ic95: bool
    conclusao_q1: str


class ResultadoMetricasModelo(BaseModel):
    """Conjunto completo de métricas obrigatórias da Seção 7.2."""

    nome_modelo: str
    estrategia_split: str
    auroc: float
    auprc: float
    sensibilidade: float
    especificidade: float
    acuracia_balanceada: float
    concordancia_categorica: float
    vme: float  # Very Major Error: predizer S quando o isolado é R
    me: float   # Major Error: predizer R quando o isolado é S
    brier_score: float
    ece: float
    cobertura_conformal: float = 1.0
    taxa_abstencao: float = 0.0
    vme_pos_abstencao: float = 0.0
    me_pos_abstencao: float = 0.0
    auprc_ic95: IntervaloConfiancaBloco | None = None
    auroc_ic95: IntervaloConfiancaBloco | None = None
    vme_ic95: IntervaloConfiancaBloco | None = None


def calcular_ece(y_true: np.ndarray, probs: np.ndarray, n_bins: int = 10) -> float:
    """Calcula o Expected Calibration Error (ECE) para avaliar a qualidade das probabilidades."""
    y_arr = np.asarray(y_true, dtype=np.float64)
    p_arr = np.asarray(probs, dtype=np.float64)
    limites = np.linspace(0.0, 1.0, n_bins + 1)
    ece = 0.0
    n = len(y_arr)
    if n == 0:
        return 0.0

    for i in range(n_bins):
        mascara = (p_arr >= limites[i]) & (p_arr < limites[i + 1] if i < n_bins - 1 else p_arr <= limites[i + 1])
        if np.any(mascara):
            acc_bin = float(np.mean(y_arr[mascara]))
            conf_bin = float(np.mean(p_arr[mascara]))
            ece += (np.sum(mascara) / n) * abs(acc_bin - conf_bin)
    return float(round(ece, 4))


def _calcular_vme_me(y_true: np.ndarray, y_pred: np.ndarray) -> tuple[float, float, float, float]:
    """
    Calcula Sensibilidade, Especificidade, VME (Very Major Error) e ME (Major Error).
    - VME = FN / (TP + FN) -> proporção de isolados verdadeiramente Resistentes classificados como Sensíveis.
    - ME  = FP / (TN + FP) -> proporção de isolados verdadeiramente Sensíveis classificados como Resistentes.
    """
    y_t = np.asarray(y_true, dtype=np.int32)
    y_p = np.asarray(y_pred, dtype=np.int32)

    mascara_r = y_t == 1
    mascara_s = y_t == 0

    tp = int(np.sum(mascara_r & (y_p == 1)))
    fn = int(np.sum(mascara_r & (y_p == 0)))
    tn = int(np.sum(mascara_s & (y_p == 0)))
    fp = int(np.sum(mascara_s & (y_p == 1)))

    total_r = tp + fn
    total_s = tn + fp

    sensibilidade = (tp / total_r) if total_r > 0 else 0.0
    especificidade = (tn / total_s) if total_s > 0 else 0.0
    vme = (fn / total_r) if total_r > 0 else 0.0
    me = (fp / total_s) if total_s > 0 else 0.0

    return float(sensibilidade), float(especificidade), float(vme), float(me)


def calcular_metricas_completas(
    nome_modelo: str,
    estrategia_split: str,
    y_true: np.ndarray,
    probs: np.ndarray,
    limiar: float = 0.50,
    predicoes_incerteza: list[PredicaoComIncerteza] | None = None,
    clusters: list[str] | np.ndarray | None = None,
    n_boot: int = 200,
    seed: int = GLOBAL_RANDOM_SEED,
) -> ResultadoMetricasModelo:
    """Calcula todas as métricas discriminativas, clínicas (VME/ME), de calibração e de abstenção."""
    y_arr = np.asarray(y_true, dtype=np.int32)
    p_arr = np.clip(np.asarray(probs, dtype=np.float64), 0.0, 1.0)
    y_pred = (p_arr >= limiar).astype(np.int32)

    if len(np.unique(y_arr)) >= 2:
        auroc = float(roc_auc_score(y_arr, p_arr))
        auprc = float(average_precision_score(y_arr, p_arr))
    else:
        auroc = 0.5
        auprc = float(np.mean(y_arr))

    sens, espec, vme, me = _calcular_vme_me(y_arr, y_pred)
    bal_acc = float(balanced_accuracy_score(y_arr, y_pred)) if len(np.unique(y_arr)) >= 2 else 0.5
    concordancia = float(np.mean(y_arr == y_pred))
    brier = float(brier_score_loss(y_arr, p_arr))
    ece = calcular_ece(y_arr, p_arr)

    cobertura = 1.0
    taxa_abst = 0.0
    vme_pos = vme
    me_pos = me

    if predicoes_incerteza is not None and len(predicoes_incerteza) == len(y_arr):
        mascara_cobertos = np.array([not p.absteve for p in predicoes_incerteza], dtype=bool)
        cobertura = float(np.mean(mascara_cobertos))
        taxa_abst = float(1.0 - cobertura)
        if np.any(mascara_cobertos):
            y_cob = y_arr[mascara_cobertos]
            pred_cob = np.array(
                [1 if p.classe == "resistente" else 0 for p in predicoes_incerteza if not p.absteve],
                dtype=np.int32,
            )
            _, _, vme_pos, me_pos = _calcular_vme_me(y_cob, pred_cob)

    auprc_ic = None
    auroc_ic = None
    vme_ic = None
    if clusters is not None and n_boot > 0:
        ics = bootstrap_em_blocos_por_cluster(
            y_true=y_arr,
            probs=p_arr,
            clusters=clusters,
            limiar=limiar,
            n_boot=n_boot,
            seed=seed,
        )
        auprc_ic = ics["auprc"]
        auroc_ic = ics["auroc"]
        vme_ic = ics["vme"]

    return ResultadoMetricasModelo(
        nome_modelo=nome_modelo,
        estrategia_split=estrategia_split,
        auroc=round(auroc, 4),
        auprc=round(auprc, 4),
        sensibilidade=round(sens, 4),
        especificidade=round(espec, 4),
        acuracia_balanceada=round(bal_acc, 4),
        concordancia_categorica=round(concordancia, 4),
        vme=round(vme, 4),
        me=round(me, 4),
        brier_score=round(brier, 4),
        ece=round(ece, 4),
        cobertura_conformal=round(cobertura, 4),
        taxa_abstencao=round(taxa_abst, 4),
        vme_pos_abstencao=round(vme_pos, 4),
        me_pos_abstencao=round(me_pos, 4),
        auprc_ic95=auprc_ic,
        auroc_ic95=auroc_ic,
        vme_ic95=vme_ic,
    )


def bootstrap_em_blocos_por_cluster(
    y_true: np.ndarray,
    probs: np.ndarray,
    clusters: list[str] | np.ndarray,
    limiar: float = 0.50,
    n_boot: int = 200,
    seed: int = GLOBAL_RANDOM_SEED,
) -> dict[str, IntervaloConfiancaBloco]:
    """
    Calcula intervalos de confiança de 95% via bootstrap em blocos (reamostrando clusters filogenéticos
    com reposição, e NÃO isolados individuais — Seção 7.2).
    """
    rng = np.random.default_rng(seed)
    y_arr = np.asarray(y_true, dtype=np.int32)
    p_arr = np.asarray(probs, dtype=np.float64)
    c_arr = np.asarray(clusters, dtype=str)

    clusters_unicos = np.unique(c_arr)
    mapa_indices = {c: np.where(c_arr == c)[0] for c in clusters_unicos}
    n_clusters = len(clusters_unicos)

    lista_auprc: list[float] = []
    lista_auroc: list[float] = []
    lista_vme: list[float] = []

    for _ in range(n_boot):
        clusters_sorteados = rng.choice(clusters_unicos, size=n_clusters, replace=True)
        idx_boot = np.concatenate([mapa_indices[c] for c in clusters_sorteados])
        yb = y_arr[idx_boot]
        pb = p_arr[idx_boot]
        if len(np.unique(yb)) < 2:
            continue
        lista_auprc.append(float(average_precision_score(yb, pb)))
        lista_auroc.append(float(roc_auc_score(yb, pb)))
        _, _, vme_b, _ = _calcular_vme_me(yb, (pb >= limiar).astype(np.int32))
        lista_vme.append(vme_b)

    def _empacotar(est: float, amostras: list[float]) -> IntervaloConfiancaBloco:
        if not amostras:
            return IntervaloConfiancaBloco(estimativa=round(est, 4), ic95_inferior=round(est, 4), ic95_superior=round(est, 4))
        inf = float(np.percentile(amostras, 2.5))
        sup = float(np.percentile(amostras, 97.5))
        return IntervaloConfiancaBloco(
            estimativa=round(est, 4),
            ic95_inferior=round(inf, 4),
            ic95_superior=round(sup, 4),
        )

    est_auprc = float(average_precision_score(y_arr, p_arr)) if len(np.unique(y_arr)) >= 2 else 0.5
    est_auroc = float(roc_auc_score(y_arr, p_arr)) if len(np.unique(y_arr)) >= 2 else 0.5
    _, _, est_vme, _ = _calcular_vme_me(y_arr, (p_arr >= limiar).astype(np.int32))

    return {
        "auprc": _empacotar(est_auprc, lista_auprc),
        "auroc": _empacotar(est_auroc, lista_auroc),
        "vme": _empacotar(est_vme, lista_vme),
    }


def comparar_modelo_vs_baseline_pareado(
    nome_modelo: str,
    y_true: np.ndarray,
    probs_modelo: np.ndarray,
    probs_b1: np.ndarray,
    clusters: list[str] | np.ndarray,
    limiar: float = 0.50,
    n_boot: int = 200,
    seed: int = GLOBAL_RANDOM_SEED,
) -> ComparacaoPareadaBaseline:
    """
    Realiza comparação pareada por bootstrap em blocos (clusters) entre um modelo de ML e o Baseline B1 (REG-04).
    O ML só é considerado superior se o limite inferior do IC 95% de Delta AUPRC (Modelo - B1) for > 0.
    """
    rng = np.random.default_rng(seed)
    y_arr = np.asarray(y_true, dtype=np.int32)
    pm_arr = np.asarray(probs_modelo, dtype=np.float64)
    pb_arr = np.asarray(probs_b1, dtype=np.float64)
    c_arr = np.asarray(clusters, dtype=str)

    clusters_unicos = np.unique(c_arr)
    mapa_indices = {c: np.where(c_arr == c)[0] for c in clusters_unicos}
    n_clusters = len(clusters_unicos)

    deltas_auprc: list[float] = []
    deltas_vme: list[float] = []

    for _ in range(n_boot):
        clusters_sorteados = rng.choice(clusters_unicos, size=n_clusters, replace=True)
        idx_boot = np.concatenate([mapa_indices[c] for c in clusters_sorteados])
        yb = y_arr[idx_boot]
        if len(np.unique(yb)) < 2:
            continue
        auprc_m = float(average_precision_score(yb, pm_arr[idx_boot]))
        auprc_b = float(average_precision_score(yb, pb_arr[idx_boot]))
        deltas_auprc.append(auprc_m - auprc_b)

        _, _, vme_m, _ = _calcular_vme_me(yb, (pm_arr[idx_boot] >= limiar).astype(np.int32))
        _, _, vme_b, _ = _calcular_vme_me(yb, (pb_arr[idx_boot] >= limiar).astype(np.int32))
        deltas_vme.append(vme_m - vme_b)

    delta_medio = float(np.mean(deltas_auprc)) if deltas_auprc else 0.0
    ic_inf = float(np.percentile(deltas_auprc, 2.5)) if deltas_auprc else 0.0
    ic_sup = float(np.percentile(deltas_auprc, 97.5)) if deltas_auprc else 0.0
    delta_vme = float(np.mean(deltas_vme)) if deltas_vme else 0.0

    supera = ic_inf > 0.0
    if supera:
        conclusao = (
            f"{nome_modelo} supera significativamente o Baseline B1 em grupos não vistos "
            f"(Delta AUPRC = +{delta_medio:.4f}, IC95% [{ic_inf:+.4f}, {ic_sup:+.4f}], Delta VME = {delta_vme:+.4f})."
        )
    else:
        conclusao = (
            f"{nome_modelo} NÃO supera com significância estatística o Baseline B1 "
            f"(Delta AUPRC = {delta_medio:+.4f}, IC95% [{ic_inf:+.4f}, {ic_sup:+.4f}]). "
            "Resultado reportado conforme REG-04."
        )

    return ComparacaoPareadaBaseline(
        nome_modelo=nome_modelo,
        delta_auprc_medio=round(delta_medio, 4),
        delta_auprc_ic95_inf=round(ic_inf, 4),
        delta_auprc_ic95_sup=round(ic_sup, 4),
        delta_vme_medio=round(delta_vme, 4),
        supera_b1_com_ic95=supera,
        conclusao_q1=conclusao,
    )

