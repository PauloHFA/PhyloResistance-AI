"""
Rastreamento de experimentos compatível com MLflow e armazenamento estruturado local (Seção 12).
Registra parâmetros de hiperparâmetros, métricas de validação por grupo, versões de bancos e artefatos.
"""

from __future__ import annotations

import json
import uuid
from dataclasses import asdict, dataclass, field
from datetime import UTC, datetime
from pathlib import Path
from typing import Any


@dataclass
class RegistroExecucaoExperimento:
    """Metadados, hiperparâmetros, métricas e artefatos de uma execução de experimento (MLflow Run)."""

    run_id: str
    experimento: str
    especie: str
    farmaco: str
    timestamp_utc: str
    parametros: dict[str, Any] = field(default_factory=dict)
    metricas: dict[str, float] = field(default_factory=dict)
    tags: dict[str, str] = field(default_factory=dict)
    artefatos: list[str] = field(default_factory=list)


class RastreadorExperimentosMLflow:
    """
    Rastreador híbrido de experimentos (Seção 12):
    - Salva estrutura compatível com MLflow local em `reports/mlruns/` (`params/`, `metrics/`, `tags/`, `meta.json`).
    - Se a biblioteca `mlflow` estiver instalada e ativa no ambiente, também sincroniza via `mlflow.log_params`/`log_metrics`.
    """

    def __init__(
        self,
        diretorio_mlruns: Path = Path("reports/mlruns"),
        nome_experimento: str = "phyloresistance-amr-benchmark",
    ) -> None:
        self.diretorio_mlruns = diretorio_mlruns
        self.nome_experimento = nome_experimento
        self.diretorio_mlruns.mkdir(parents=True, exist_ok=True)

    def registrar_execucao_farmaco(
        self,
        especie: str,
        farmaco: str,
        parametros: dict[str, Any],
        metricas: dict[str, float],
        tags: dict[str, str] | None = None,
        artefatos: list[str] | None = None,
    ) -> RegistroExecucaoExperimento:
        """Cria um registro completo da execução de treinamento/validação de um par bactéria-fármaco."""
        run_id = uuid.uuid4().hex[:12]
        agora = datetime.now(UTC).isoformat()
        tags_finais = {
            "especie": especie,
            "farmaco": farmaco,
            "framework": "PhyloResistance-AI",
            **(tags or {}),
        }
        registro = RegistroExecucaoExperimento(
            run_id=run_id,
            experimento=self.nome_experimento,
            especie=especie,
            farmaco=farmaco,
            timestamp_utc=agora,
            parametros=dict(parametros),
            metricas={k: round(float(v), 4) for k, v in metricas.items()},
            tags=tags_finais,
            artefatos=list(artefatos or []),
        )

        # Estrutura compatível com MLflow FileStore
        pasta_exp = self.diretorio_mlruns / self.nome_experimento
        pasta_run = pasta_exp / run_id
        (pasta_run / "params").mkdir(parents=True, exist_ok=True)
        (pasta_run / "metrics").mkdir(parents=True, exist_ok=True)
        (pasta_run / "tags").mkdir(parents=True, exist_ok=True)

        for k, v in registro.parametros.items():
            (pasta_run / "params" / k).write_text(str(v), encoding="utf-8")
        for k, v in registro.metricas.items():
            (pasta_run / "metrics" / k).write_text(f"{agora} {v}\n", encoding="utf-8")
        for k, v in registro.tags.items():
            (pasta_run / "tags" / k).write_text(str(v), encoding="utf-8")

        (pasta_run / "meta.json").write_text(
            json.dumps(asdict(registro), ensure_ascii=False, indent=2),
            encoding="utf-8",
        )

        # Atualiza índice consolidado de runs
        caminho_indice = self.diretorio_mlruns / "indice_experimentos.json"
        historico: list[dict[str, Any]] = []
        if caminho_indice.exists():
            try:
                historico = json.loads(caminho_indice.read_text(encoding="utf-8"))
            except Exception:
                historico = []
        historico.append(asdict(registro))
        caminho_indice.write_text(
            json.dumps(historico, ensure_ascii=False, indent=2),
            encoding="utf-8",
        )

        # Integração opcional com pacote mlflow se instalado
        try:
            import mlflow  # type: ignore[import-not-found]

            mlflow.set_tracking_uri(self.diretorio_mlruns.resolve().as_uri())
            mlflow.set_experiment(self.nome_experimento)
            with mlflow.start_run(run_name=f"{especie}_{farmaco}"):
                mlflow.log_params({k: str(v) for k, v in registro.parametros.items()})
                mlflow.log_metrics(registro.metricas)
                mlflow.set_tags(registro.tags)
        except Exception:
            pass

        return registro
