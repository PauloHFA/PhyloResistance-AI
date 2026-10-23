"""Escada de modelos B0 a B5, controle negativo de k-mers curtos e ajuste Optuna por grupo (Seção 6)."""

from __future__ import annotations

import numpy as np
import optuna
import scipy.sparse as sp
import xgboost as xgb
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score
from sklearn.model_selection import StratifiedGroupKFold

from phyloresistance.config import (
    CATALOGO_DETERMINANTES_ECOLI,
    GLOBAL_RANDOM_SEED,
)
from phyloresistance.genomica.kmers import FiltroSelecaoKmersDobra

optuna.logging.set_verbosity(optuna.logging.WARNING)


def calcular_peso_classe_balanceado(y_treino: np.ndarray) -> float:
    """Calcula scale_pos_weight = N_negativos / N_positivos apenas nos dados de treino (Seção 6.2)."""
    y_arr = np.asarray(y_treino, dtype=np.int32)
    n_pos = int(np.sum(y_arr == 1))
    n_neg = int(np.sum(y_arr == 0))
    if n_pos == 0:
        return 1.0
    return float(max(0.1, n_neg / n_pos))


class ClassificadorMajoritarioB0:
    """Nível B0: Classificador trivial de classe majoritária (piso de comparação)."""

    def __init__(self) -> None:
        self.prob_positiva_: float = 0.5
        self.classe_majoritaria_: int = 0

    def fit(self, X: np.ndarray | sp.csr_matrix, y: np.ndarray) -> ClassificadorMajoritarioB0:
        y_arr = np.asarray(y, dtype=np.int32)
        self.prob_positiva_ = float(np.mean(y_arr)) if len(y_arr) > 0 else 0.5
        self.classe_majoritaria_ = 1 if self.prob_positiva_ >= 0.5 else 0
        return self

    def predict_proba(self, X: np.ndarray | sp.csr_matrix) -> np.ndarray:
        n = X.shape[0]
        p1 = np.full(n, self.prob_positiva_, dtype=np.float64)
        return np.column_stack([1.0 - p1, p1])

    def predict(self, X: np.ndarray | sp.csr_matrix) -> np.ndarray:
        return np.full(X.shape[0], self.classe_majoritaria_, dtype=np.int32)


class BaselineRegrasB1:
    """
    Nível B1: Baseline Científico Obrigatório por Regras (CARD / AMRFinderPlus — Seção 6.1 e REG-04).
    Prediz Resistente (1) se houver pelo menos 1 determinante conhecido (F1 ou F2)
    catalogado para o fármaco avaliado.
    """

    def __init__(self, farmaco: str, nomes_atributos_f1_f2: list[str]) -> None:
        self.farmaco = farmaco
        self.nomes_atributos_f1_f2 = list(nomes_atributos_f1_f2)
        self.indices_determinantes_farmaco_: list[int] = []
        self.taxa_base_sem_det_: float = 0.05
        self.taxa_com_det_: float = 0.92

    def fit(self, X_f1_f2: np.ndarray, y: np.ndarray) -> BaselineRegrasB1:
        determinantes_validos = {
            d.id_atributo
            for d in CATALOGO_DETERMINANTES_ECOLI
            if self.farmaco in d.farmacos_associados
        }
        self.indices_determinantes_farmaco_ = [
            idx
            for idx, nome in enumerate(self.nomes_atributos_f1_f2)
            if nome in determinantes_validos
        ]
        X_arr = np.asarray(X_f1_f2, dtype=np.float64)
        y_arr = np.asarray(y, dtype=np.int32)

        if self.indices_determinantes_farmaco_:
            contagem = X_arr[:, self.indices_determinantes_farmaco_].sum(axis=1)
            tem_det = contagem > 0
            if np.any(~tem_det):
                self.taxa_base_sem_det_ = float(np.clip(np.mean(y_arr[~tem_det]), 0.02, 0.45))
            if np.any(tem_det):
                self.taxa_com_det_ = float(np.clip(np.mean(y_arr[tem_det]), 0.55, 0.98))
        return self

    def predict(self, X_f1_f2: np.ndarray) -> np.ndarray:
        X_arr = np.asarray(X_f1_f2, dtype=np.float64)
        if not self.indices_determinantes_farmaco_:
            return np.zeros(X_arr.shape[0], dtype=np.int32)
        contagem = X_arr[:, self.indices_determinantes_farmaco_].sum(axis=1)
        return (contagem >= 1.0).astype(np.int32)

    def predict_proba(self, X_f1_f2: np.ndarray) -> np.ndarray:
        X_arr = np.asarray(X_f1_f2, dtype=np.float64)
        if not self.indices_determinantes_farmaco_:
            p1 = np.full(X_arr.shape[0], self.taxa_base_sem_det_)
            return np.column_stack([1.0 - p1, p1])

        contagem = X_arr[:, self.indices_determinantes_farmaco_].sum(axis=1)
        # Escala monotônica conforme o número de determinantes conhecidos presentes
        p1 = np.where(
            contagem >= 1.0,
            np.clip(self.taxa_com_det_ + 0.03 * (contagem - 1.0), 0.51, 0.99),
            self.taxa_base_sem_det_,
        )
        return np.column_stack([1.0 - p1, p1])


class ModeloLogisticoRegularizadoB2:
    """Nível B2: Regressão logística regularizada L1 (Lasso) sobre F1 + F2."""

    def __init__(self, C: float = 1.0, seed: int = GLOBAL_RANDOM_SEED) -> None:
        self.C = C
        self.seed = seed
        self.modelo_ = LogisticRegression(
            l1_ratio=1.0,
            solver="liblinear",
            C=C,
            class_weight="balanced",
            random_state=seed,
            max_iter=1000,
        )

    def fit(self, X_f1_f2: np.ndarray, y: np.ndarray) -> ModeloLogisticoRegularizadoB2:
        self.modelo_.fit(np.asarray(X_f1_f2, dtype=np.float64), np.asarray(y, dtype=np.int32))
        return self

    def predict_proba(self, X_f1_f2: np.ndarray) -> np.ndarray:
        return self.modelo_.predict_proba(np.asarray(X_f1_f2, dtype=np.float64))

    def predict(self, X_f1_f2: np.ndarray) -> np.ndarray:
        return (self.predict_proba(X_f1_f2)[:, 1] >= 0.5).astype(np.int32)


class ModeloArvoreCuradoB3:
    """Nível B3: XGBoost sobre F1 + F2 (determinantes curados)."""

    def __init__(
        self,
        n_estimators: int = 200,
        max_depth: int = 4,
        learning_rate: float = 0.05,
        seed: int = GLOBAL_RANDOM_SEED,
    ) -> None:
        self.n_estimators = n_estimators
        self.max_depth = max_depth
        self.learning_rate = learning_rate
        self.seed = seed
        self.modelo_: xgb.XGBClassifier | None = None

    def fit(self, X_f1_f2: np.ndarray, y: np.ndarray) -> ModeloArvoreCuradoB3:
        y_arr = np.asarray(y, dtype=np.int32)
        peso = calcular_peso_classe_balanceado(y_arr)
        self.modelo_ = xgb.XGBClassifier(
            n_estimators=self.n_estimators,
            max_depth=self.max_depth,
            learning_rate=self.learning_rate,
            subsample=0.8,
            colsample_bytree=0.8,
            scale_pos_weight=peso,
            eval_metric="aucpr",
            random_state=self.seed,
            n_jobs=1,
        )
        self.modelo_.fit(np.asarray(X_f1_f2, dtype=np.float64), y_arr)
        return self

    def predict_proba(self, X_f1_f2: np.ndarray) -> np.ndarray:
        if self.modelo_ is None:
            raise RuntimeError("Modelo B3 não ajustado.")
        return self.modelo_.predict_proba(np.asarray(X_f1_f2, dtype=np.float64))

    def predict(self, X_f1_f2: np.ndarray) -> np.ndarray:
        return (self.predict_proba(X_f1_f2)[:, 1] >= 0.5).astype(np.int32)


class PipelineCompletoB4:
    """
    Nível B4: Modelo Completo XGBoost sobre F1 + F2 + F3 (captura variantes novas e interações).
    Cumpre estritamente a Regra REG-03: o filtro de frequência de k-mers e a seleção supervisionada
    em F3 são ajustados APENAS nos dados de treino desta dobra.
    """

    def __init__(
        self,
        nomes_f1_f2: list[str],
        nomes_f3: list[str],
        max_kmers_f3: int = 60,
        min_freq_kmer: float = 0.02,
        max_freq_kmer: float = 0.98,
        n_estimators: int = 250,
        max_depth: int = 4,
        learning_rate: float = 0.05,
        seed: int = GLOBAL_RANDOM_SEED,
    ) -> None:
        self.nomes_f1_f2 = list(nomes_f1_f2)
        self.nomes_f3 = list(nomes_f3)
        self.filtro_f3 = FiltroSelecaoKmersDobra(
            min_freq=min_freq_kmer,
            max_freq=max_freq_kmer,
            max_kmers_selecionados=max_kmers_f3,
        )
        self.n_estimators = n_estimators
        self.max_depth = max_depth
        self.learning_rate = learning_rate
        self.seed = seed
        self.modelo_: xgb.XGBClassifier | None = None
        self.nomes_atributos_finais_: list[str] = []

    def _combinar_matrizes(
        self, X_f1_f2: np.ndarray, X_f3_filtrado: sp.csr_matrix
    ) -> np.ndarray:
        X_curado = np.asarray(X_f1_f2, dtype=np.float64)
        X_kmers = X_f3_filtrado.toarray()
        return np.hstack([X_curado, X_kmers])

    def fit(
        self,
        X_f1_f2_treino: np.ndarray,
        X_f3_treino: sp.csr_matrix | np.ndarray,
        y_treino: np.ndarray,
        ids_treino: list[str] | None = None,
    ) -> PipelineCompletoB4:
        y_arr = np.asarray(y_treino, dtype=np.int32)
        X_f3_sel = self.filtro_f3.fit_transform(X_f3_treino, y_arr, ids_treino=ids_treino)
        nomes_f3_sel = self.filtro_f3.filtrar_nomes(self.nomes_f3)
        self.nomes_atributos_finais_ = self.nomes_f1_f2 + nomes_f3_sel

        X_treino_full = self._combinar_matrizes(X_f1_f2_treino, X_f3_sel)
        peso = calcular_peso_classe_balanceado(y_arr)

        self.modelo_ = xgb.XGBClassifier(
            n_estimators=self.n_estimators,
            max_depth=self.max_depth,
            learning_rate=self.learning_rate,
            subsample=0.8,
            colsample_bytree=0.6,
            scale_pos_weight=peso,
            eval_metric="aucpr",
            random_state=self.seed,
            n_jobs=1,
        )
        self.modelo_.fit(X_treino_full, y_arr)
        return self

    def transformar_atributos(
        self,
        X_f1_f2: np.ndarray,
        X_f3: sp.csr_matrix | np.ndarray,
    ) -> np.ndarray:
        X_f3_sel = self.filtro_f3.transform(X_f3)
        return self._combinar_matrizes(X_f1_f2, X_f3_sel)

    def predict_proba(
        self,
        X_f1_f2: np.ndarray,
        X_f3: sp.csr_matrix | np.ndarray,
    ) -> np.ndarray:
        if self.modelo_ is None:
            raise RuntimeError("PipelineCompletoB4 não ajustado.")
        X_full = self.transformar_atributos(X_f1_f2, X_f3)
        return self.modelo_.predict_proba(X_full)

    def predict(
        self,
        X_f1_f2: np.ndarray,
        X_f3: sp.csr_matrix | np.ndarray,
    ) -> np.ndarray:
        return (self.predict_proba(X_f1_f2, X_f3)[:, 1] >= 0.5).astype(np.int32)


class ModeloLinearEsparsoKmersB5:
    """Nível B5: Regressão logística regularizada L1 apenas sobre F3 (k-mers longos filtrados na dobra)."""

    def __init__(
        self,
        max_kmers_f3: int = 80,
        C: float = 0.8,
        seed: int = GLOBAL_RANDOM_SEED,
    ) -> None:
        self.filtro_f3 = FiltroSelecaoKmersDobra(max_kmers_selecionados=max_kmers_f3)
        self.modelo_ = LogisticRegression(
            l1_ratio=1.0,
            solver="liblinear",
            C=C,
            class_weight="balanced",
            random_state=seed,
            max_iter=1000,
        )

    def fit(
        self,
        X_f3_treino: sp.csr_matrix | np.ndarray,
        y_treino: np.ndarray,
        ids_treino: list[str] | None = None,
    ) -> ModeloLinearEsparsoKmersB5:
        y_arr = np.asarray(y_treino, dtype=np.int32)
        X_sel = self.filtro_f3.fit_transform(X_f3_treino, y_arr, ids_treino=ids_treino)
        self.modelo_.fit(X_sel, y_arr)
        return self

    def predict_proba(self, X_f3: sp.csr_matrix | np.ndarray) -> np.ndarray:
        X_sel = self.filtro_f3.transform(X_f3)
        return self.modelo_.predict_proba(X_sel)

    def predict(self, X_f3: sp.csr_matrix | np.ndarray) -> np.ndarray:
        return (self.predict_proba(X_f3)[:, 1] >= 0.5).astype(np.int32)


class ModeloControleKmersCurtos:
    """
    Controle Negativo (Ablação Seção 2 e 6.1):
    XGBoost treinado sobre frequências de k-mers curtos (k=6 a 8).
    Demonstra que frequência de k-mer curto mede apenas composição genômica global, não resistência.
    """

    def __init__(self, seed: int = GLOBAL_RANDOM_SEED) -> None:
        self.seed = seed
        self.modelo_: xgb.XGBClassifier | None = None

    def fit(self, X_k_curto: np.ndarray, y: np.ndarray) -> ModeloControleKmersCurtos:
        y_arr = np.asarray(y, dtype=np.int32)
        peso = calcular_peso_classe_balanceado(y_arr)
        self.modelo_ = xgb.XGBClassifier(
            n_estimators=100,
            max_depth=3,
            learning_rate=0.05,
            scale_pos_weight=peso,
            eval_metric="aucpr",
            random_state=self.seed,
            n_jobs=1,
        )
        self.modelo_.fit(np.asarray(X_k_curto, dtype=np.float64), y_arr)
        return self

    def predict_proba(self, X_k_curto: np.ndarray) -> np.ndarray:
        if self.modelo_ is None:
            raise RuntimeError("ModeloControleKmersCurtos não ajustado.")
        return self.modelo_.predict_proba(np.asarray(X_k_curto, dtype=np.float64))

    def predict(self, X_k_curto: np.ndarray) -> np.ndarray:
        return (self.predict_proba(X_k_curto)[:, 1] >= 0.5).astype(np.int32)


def otimizar_xgb_optuna_por_grupo(
    X_treino: np.ndarray,
    y_treino: np.ndarray,
    grupos_treino: np.ndarray | list[str],
    n_trials: int = 8,
    seed: int = GLOBAL_RANDOM_SEED,
) -> dict[str, int | float]:
    """
    Ajuste de hiperparâmetros com Optuna em validação aninhada por grupo filogenético (Seção 6.2 e REG-03).
    Usa APENAS os dados e grupos da dobra de treino.
    """
    X_arr = np.asarray(X_treino, dtype=np.float64)
    y_arr = np.asarray(y_treino, dtype=np.int32)
    g_arr = np.asarray(grupos_treino, dtype=str)

    grupos_unicos = np.unique(g_arr)
    if len(grupos_unicos) < 3:
        return {"n_estimators": 200, "max_depth": 4, "learning_rate": 0.05}

    cv_interno = StratifiedGroupKFold(n_splits=3, shuffle=True, random_state=seed)
    peso = calcular_peso_classe_balanceado(y_arr)

    def objective(trial: optuna.Trial) -> float:
        n_estimators = trial.suggest_int("n_estimators", 120, 300, step=60)
        max_depth = trial.suggest_int("max_depth", 3, 6)
        learning_rate = trial.suggest_float("learning_rate", 0.02, 0.12, log=True)

        scores: list[float] = []
        for idx_tr, idx_val in cv_interno.split(X_arr, y_arr, groups=g_arr):
            if len(np.unique(y_arr[idx_val])) < 2:
                continue
            clf = xgb.XGBClassifier(
                n_estimators=n_estimators,
                max_depth=max_depth,
                learning_rate=learning_rate,
                subsample=0.8,
                colsample_bytree=0.7,
                scale_pos_weight=peso,
                eval_metric="aucpr",
                random_state=seed,
                n_jobs=1,
            )
            clf.fit(X_arr[idx_tr], y_arr[idx_tr])
            prob = clf.predict_proba(X_arr[idx_val])[:, 1]
            scores.append(float(average_precision_score(y_arr[idx_val], prob)))
        return float(np.mean(scores)) if scores else 0.5

    sampler = optuna.samplers.TPESampler(seed=seed)
    estudo = optuna.create_study(direction="maximize", sampler=sampler)
    estudo.optimize(objective, n_trials=n_trials, show_progress_bar=False)
    return dict(estudo.best_params)
