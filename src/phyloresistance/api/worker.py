"""
Worker assíncrono para processamento pesado de genomas (QC, k-mers, predição e SHAP)
com adaptador para Redis Queue (RQ) e fallback em ThreadPool local (Seção 9.1 e Fase 5/7).
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any

from phyloresistance.api.app import _EXECUTOR_FILA, processar_genoma_worker


class GerenciadorFilaWorkers:
    """
    Abstrai o despacho de jobs de predição genômica:
    - Se REDIS_URL estiver configurada e biblioteca `rq`/`redis` disponível, enfileira no Redis Queue.
    - Caso contrário, despacha no ThreadPoolExecutor local sem bloqueio de dependência externa.
    """

    def __init__(self, redis_url: str | None = None, nome_fila: str = "phyloresistance_jobs") -> None:
        self.redis_url = redis_url or os.environ.get("REDIS_URL")
        self.nome_fila = nome_fila
        self._fila_rq: Any | None = None
        if self.redis_url:
            try:
                import redis  # type: ignore[import-not-found]
                from rq import Queue  # type: ignore[import-not-found]

                conexao = redis.from_url(self.redis_url)
                self._fila_rq = Queue(self.nome_fila, connection=conexao)
            except Exception:
                self._fila_rq = None

    @property
    def modo_execucao(self) -> str:
        return "rq_redis" if self._fila_rq is not None else "threadpool_local"

    def enfileirar_job_genoma(
        self,
        job_id: str,
        caminho_entrada: Path,
        formato_entrada: str = "FASTA",
        especie: str = "Escherichia coli",
        farmacos_solicitados: list[str] | None = None,
    ) -> dict[str, Any]:
        """Despacha o processamento pesado do isolado para o worker e retorna metadados de fila."""
        if self._fila_rq is not None:
            job_rq = self._fila_rq.enqueue(
                processar_genoma_worker,
                job_id,
                caminho_entrada,
                formato_entrada,
                especie,
                farmacos_solicitados,
                job_id=job_id,
                job_timeout=180,
            )
            return {"job_id": job_id, "modo": "rq_redis", "rq_id": getattr(job_rq, "id", job_id)}

        futuro = _EXECUTOR_FILA.submit(
            processar_genoma_worker,
            job_id,
            caminho_entrada,
            formato_entrada,
            especie,
            farmacos_solicitados,
        )
        try:
            futuro.result(timeout=3.0)
        except Exception:
            pass
        return {"job_id": job_id, "modo": "threadpool_local"}

