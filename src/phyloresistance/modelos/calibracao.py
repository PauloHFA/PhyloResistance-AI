"""Calibração probabilística (Isotônica / Platt) e Predição Conformal com Abstenção (Seção 6.2, Q4, REG-11)."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np
from sklearn.isotonic import IsotonicRegression
from sklearn.linear_model import LogisticRegression


@dataclass(frozen=True)
class PredicaoComIncerteza:
    """Resultado individual de predição calibrada com conjunto conformal e abstenção."""

    probabilidade_resistencia: float
    classe: Literal["resistente", "sensivel", "indeterminado"]
    confianca: Literal["alta", "media", "baixa"]
    conjunto_conformal: tuple[int, ...]
    absteve: bool


class CalibradorConformalAMR:
    """
    Realiza calibração probabilística (Platt ou Isotônica) em conjunto de calibração separado
    e aplica predição conformal indutiva (tipo MAPIE / split conformal) para devolver
    'indeterminado' quando a incerteza for alta, reduzindo erros graves (VME e ME).
    """

    def __init__(
        self,
        metodo: Literal["platt", "isotonica"] = "platt",
        alpha: float = 0.12,
        limiar_decisao: float = 0.50,
    ) -> None:
        self.metodo = metodo
        self.alpha = alpha
        self.limiar_decisao = limiar_decisao
        self._platt: LogisticRegression | None = None
        self._iso: IsotonicRegression | None = None
        self.quantil_nao_conformidade_: float = 0.38
        self.faixa_abstencao_: tuple[float, float] = (0.35, 0.65)

    def fit(self, probs_brutas_cal: np.ndarray, y_cal: np.ndarray) -> CalibradorConformalAMR:
        """Ajusta o calibrador e os limiares conformais em dados de calibração separados."""
        p_arr = np.clip(np.asarray(probs_brutas_cal, dtype=np.float64).ravel(), 1e-5, 1.0 - 1e-5)
        y_arr = np.asarray(y_cal, dtype=np.int32).ravel()

        if len(np.unique(y_arr)) >= 2 and len(y_arr) >= 10:
            if self.metodo == "isotonica":
                self._iso = IsotonicRegression(out_of_bounds="clip", y_min=0.01, y_max=0.99)
                self._iso.fit(p_arr, y_arr)
            else:
                logits = np.log(p_arr / (1.0 - p_arr)).reshape(-1, 1)
                self._platt = LogisticRegression(C=1.0, random_state=42)
                self._platt.fit(logits, y_arr)

        probs_calibradas = self.calibrar_probabilidades(p_arr)

        # Escores de não-conformidade para a classe verdadeira: s_i = 1 - p_cal(y_i)
        prob_classe_verdadeira = np.where(y_arr == 1, probs_calibradas, 1.0 - probs_calibradas)
        escores_nc = 1.0 - prob_classe_verdadeira
        n = len(escores_nc)
        nivel_q = min(0.99, np.ceil((n + 1) * (1.0 - self.alpha)) / max(1, n))
        q_hat = float(np.quantile(escores_nc, nivel_q))
        # Define faixa de abstenção simétrica ao redor do limiar caso o conjunto conformal seja ambíguo
        self.quantil_nao_conformidade_ = float(np.clip(q_hat, 0.25, 0.45))
        margem = max(0.12, 0.50 - (1.0 - self.quantil_nao_conformidade_))
        margem = min(0.22, margem)
        self.faixa_abstencao_ = (
            round(self.limiar_decisao - margem, 4),
            round(self.limiar_decisao + margem, 4),
        )
        return self

    def calibrar_probabilidades(self, probs_brutas: np.ndarray) -> np.ndarray:
        """Transforma probabilidades brutas em probabilidades calibradas."""
        p_arr = np.clip(np.asarray(probs_brutas, dtype=np.float64).ravel(), 1e-5, 1.0 - 1e-5)
        if self._iso is not None:
            return np.clip(self._iso.predict(p_arr), 0.01, 0.99)
        if self._platt is not None:
            logits = np.log(p_arr / (1.0 - p_arr)).reshape(-1, 1)
            return np.clip(self._platt.predict_proba(logits)[:, 1], 0.01, 0.99)
        return p_arr

    def predizer_com_incerteza(
        self, probs_brutas: np.ndarray
    ) -> list[PredicaoComIncerteza]:
        """
        Gera predições calibradas com conjunto conformal e abstenção explícita ('indeterminado')
        para isolados na zona de dúvida (Seção 6.2 e Q4).
        """
        probs_cal = self.calibrar_probabilidades(probs_brutas)
        lim_inf, lim_sup = self.faixa_abstencao_
        resultados: list[PredicaoComIncerteza] = []

        for p in probs_cal:
            p_float = float(round(p, 4))
            # Verifica se cai na zona de incerteza conformal
            if lim_inf <= p_float <= lim_sup:
                resultados.append(
                    PredicaoComIncerteza(
                        probabilidade_resistencia=p_float,
                        classe="indeterminado",
                        confianca="baixa",
                        conjunto_conformal=(0, 1),
                        absteve=True,
                    )
                )
            elif p_float > lim_sup:
                conf: Literal["alta", "media"] = "alta" if p_float >= 0.80 else "media"
                resultados.append(
                    PredicaoComIncerteza(
                        probabilidade_resistencia=p_float,
                        classe="resistente",
                        confianca=conf,
                        conjunto_conformal=(1,),
                        absteve=False,
                    )
                )
            else:
                conf = "alta" if p_float <= 0.20 else "media"
                resultados.append(
                    PredicaoComIncerteza(
                        probabilidade_resistencia=p_float,
                        classe="sensivel",
                        confianca=conf,
                        conjunto_conformal=(0,),
                        absteve=False,
                    )
                )
        return resultados

