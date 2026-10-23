"""Extensões de modelagem da Fase 6: Regressão log2(MIC), Modelo Neural de Sequência B6 e Modelo Multi-Saída (Seção 2, 6.1, 6.2)."""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp
import xgboost as xgb
from sklearn.multioutput import MultiOutputClassifier
from sklearn.neural_network import MLPClassifier

from phyloresistance.config import GLOBAL_RANDOM_SEED
from phyloresistance.genomica.kmers import FiltroSelecaoKmersDobra


class RegressorMICLog2:
    """
    Regressão quantitativa de log2(MIC em mg/L) sobre F1 + F2 + F3 (Seção 6.2 e Fase 6).
    Permite avaliar Concordância Essencial (erro de até +-1 diluição dupla segundo ISO 20776-2).
    """

    def __init__(
        self,
        max_kmers_f3: int = 50,
        n_estimators: int = 200,
        max_depth: int = 4,
        learning_rate: float = 0.05,
        seed: int = GLOBAL_RANDOM_SEED,
    ) -> None:
        self.filtro_f3 = FiltroSelecaoKmersDobra(max_kmers_selecionados=max_kmers_f3)
        self.n_estimators = n_estimators
        self.max_depth = max_depth
        self.learning_rate = learning_rate
        self.seed = seed
        self.regressor_: xgb.XGBRegressor | None = None

    def fit(
        self,
        X_f1_f2_treino: np.ndarray,
        X_f3_treino: sp.csr_matrix | np.ndarray,
        mic_mg_l_treino: np.ndarray,
    ) -> RegressorMICLog2:
        mic_arr = np.clip(np.asarray(mic_mg_l_treino, dtype=np.float64), 1e-4, 1024.0)
        y_log2 = np.log2(mic_arr)
        y_bin_proxy = (y_log2 >= np.median(y_log2)).astype(np.int32)

        X_f3_sel = self.filtro_f3.fit_transform(X_f3_treino, y_bin_proxy)
        X_full = np.hstack([np.asarray(X_f1_f2_treino, dtype=np.float64), X_f3_sel.toarray()])

        self.regressor_ = xgb.XGBRegressor(
            n_estimators=self.n_estimators,
            max_depth=self.max_depth,
            learning_rate=self.learning_rate,
            subsample=0.8,
            colsample_bytree=0.7,
            random_state=self.seed,
            n_jobs=1,
        )
        self.regressor_.fit(X_full, y_log2)
        return self

    def predict_log2_mic(
        self,
        X_f1_f2: np.ndarray,
        X_f3: sp.csr_matrix | np.ndarray,
    ) -> np.ndarray:
        if self.regressor_ is None:
            raise RuntimeError("RegressorMICLog2 não ajustado.")
        X_f3_sel = self.filtro_f3.transform(X_f3)
        X_full = np.hstack([np.asarray(X_f1_f2, dtype=np.float64), X_f3_sel.toarray()])
        return np.asarray(self.regressor_.predict(X_full), dtype=np.float64)

    def predict_mic_mg_l(
        self,
        X_f1_f2: np.ndarray,
        X_f3: sp.csr_matrix | np.ndarray,
        arredondar_diluicao_dupla: bool = True,
    ) -> np.ndarray:
        log2_pred = self.predict_log2_mic(X_f1_f2, X_f3)
        if arredondar_diluicao_dupla:
            return np.power(2.0, np.round(log2_pred))
        return np.power(2.0, log2_pred)


class ModeloSequenciaConvolucionalB6:
    """
    Nível B6 da Escada de Modelos (Seção 6.1 e Fase 6):
    Extrator convolucional 1D (banco de filtros de janela deslizante sobre loci/sequência)
    seguido de Rede Neural Multicamadas (MLP) com regularização L2 e early-stopping.
    """

    def __init__(
        self,
        n_filtros_conv: int = 24,
        tamanho_kernel: int = 5,
        hidden_layer_sizes: tuple[int, ...] = (48, 24),
        seed: int = GLOBAL_RANDOM_SEED,
    ) -> None:
        self.n_filtros_conv = n_filtros_conv
        self.tamanho_kernel = tamanho_kernel
        self.hidden_layer_sizes = hidden_layer_sizes
        self.seed = seed
        self.filtro_f3 = FiltroSelecaoKmersDobra(max_kmers_selecionados=48)
        self.kernels_conv_: np.ndarray | None = None
        self.mlp_: MLPClassifier | None = None

    def _aplicar_banco_convolucional_1d(self, X_denso: np.ndarray) -> np.ndarray:
        """Aplica convolução 1D com ativação ReLU e Global Max/Mean Pooling combinada aos sinais diretos."""
        assert self.kernels_conv_ is not None
        n_amostras, largura = X_denso.shape
        k = min(self.tamanho_kernel, largura)
        n_pos = max(1, largura - k + 1)

        # Extrai janelas deslizantes (n_amostras, n_pos, k)
        janelas = np.stack([X_denso[:, i : i + k] for i in range(n_pos)], axis=1)
        # Projeção convolucional (n_amostras, n_pos, n_filtros)
        ativacoes = np.maximum(0.0, np.einsum("npk,fk->npf", janelas, self.kernels_conv_[:, :k]))
        max_pool = np.max(ativacoes, axis=1)
        mean_pool = np.mean(ativacoes, axis=1)
        return np.hstack([X_denso, max_pool, mean_pool])

    def fit(
        self,
        X_f1_f2_treino: np.ndarray,
        X_f3_treino: sp.csr_matrix | np.ndarray,
        y_treino: np.ndarray,
    ) -> ModeloSequenciaConvolucionalB6:
        rng = np.random.default_rng(self.seed)
        y_arr = np.asarray(y_treino, dtype=np.int32)
        X_f3_sel = self.filtro_f3.fit_transform(X_f3_treino, y_arr)
        X_base = np.hstack([np.asarray(X_f1_f2_treino, dtype=np.float64), X_f3_sel.toarray()])

        k = min(self.tamanho_kernel, X_base.shape[1])
        self.kernels_conv_ = rng.normal(0.0, 0.35, size=(self.n_filtros_conv, k))
        X_conv = self._aplicar_banco_convolucional_1d(X_base)

        self.mlp_ = MLPClassifier(
            hidden_layer_sizes=self.hidden_layer_sizes,
            activation="relu",
            alpha=0.05,
            max_iter=350,
            random_state=self.seed,
        )
        self.mlp_.fit(X_conv, y_arr)
        return self

    def predict_proba(
        self,
        X_f1_f2: np.ndarray,
        X_f3: sp.csr_matrix | np.ndarray,
    ) -> np.ndarray:
        if self.mlp_ is None:
            raise RuntimeError("ModeloSequenciaConvolucionalB6 não ajustado.")
        X_f3_sel = self.filtro_f3.transform(X_f3)
        X_base = np.hstack([np.asarray(X_f1_f2, dtype=np.float64), X_f3_sel.toarray()])
        X_conv = self._aplicar_banco_convolucional_1d(X_base)
        return self.mlp_.predict_proba(X_conv)

    def predict(
        self,
        X_f1_f2: np.ndarray,
        X_f3: sp.csr_matrix | np.ndarray,
    ) -> np.ndarray:
        return (self.predict_proba(X_f1_f2, X_f3)[:, 1] >= 0.5).astype(np.int32)


class ModeloMultiSaidaAMR:
    """
    Extensão Multi-Saída Opcional (Seção 2):
    Prediz simultaneamente o vetor de resistência para múltiplos fármacos da mesma espécie.
    """

    def __init__(self, seed: int = GLOBAL_RANDOM_SEED) -> None:
        self.seed = seed
        base_xgb = xgb.XGBClassifier(
            n_estimators=150,
            max_depth=4,
            learning_rate=0.05,
            subsample=0.8,
            eval_metric="logloss",
            random_state=seed,
            n_jobs=1,
        )
        self.modelo_ = MultiOutputClassifier(base_xgb)

    def fit(self, X: np.ndarray, Y_multifarmaco: np.ndarray) -> ModeloMultiSaidaAMR:
        self.modelo_.fit(
            np.asarray(X, dtype=np.float64),
            np.asarray(Y_multifarmaco, dtype=np.int32),
        )
        return self

    def predict(self, X: np.ndarray) -> np.ndarray:
        return np.asarray(self.modelo_.predict(np.asarray(X, dtype=np.float64)), dtype=np.int32)

    def predict_proba_matriz(self, X: np.ndarray) -> np.ndarray:
        """Retorna matriz (n_amostras, n_farmacos) com as probabilidades de resistência."""
        lista_probs = self.modelo_.predict_proba(np.asarray(X, dtype=np.float64))
        colunas = [p[:, 1] for p in lista_probs]
        return np.column_stack(colunas)

