"""Esquemas Pydantic de entrada, saída, MIC e observabilidade/drift da API FastAPI (Seção 9.2, 9.3, 9.4, REG-14, REG-15)."""

from __future__ import annotations

from typing import Literal

from pydantic import BaseModel, Field

from phyloresistance.config import AVISO_USO_PESQUISA


class RespostaJobCriado(BaseModel):
    """Resposta imediata (HTTP 202 Accepted) ao enfileirar um genoma para predição."""

    job_id: str
    status: Literal["em_fila", "processando", "concluido"] = "em_fila"
    formato_entrada: Literal["FASTA", "FASTQ"] = "FASTA"


class ResumoQCResposta(BaseModel):
    """Resumo das métricas de QC do genoma avaliado no serviço."""

    tamanho_pb: int
    contigs: int
    n50: int
    aprovado: bool
    motivos_reprovacao: list[str] = Field(default_factory=list)


class ItemEvidenciaAPI(BaseModel):
    """Evidência genômica individual retornada na resposta da API (Seção 9.3)."""

    tipo: str
    gene: str
    alteracao: str | None = None
    fonte: str
    contribuicao_shap: float
    status: Literal["conhecida", "candidata"]


class ItemPredicaoFarmacoAPI(BaseModel):
    """Predição por fármaco para um isolado, incluindo estimativa quantitativa opcional de MIC."""

    farmaco: str
    probabilidade_resistencia: float | None = None
    classe: Literal["resistente", "sensivel", "indeterminado", "nao_suportado"]
    confianca: Literal["alta", "media", "baixa", "n/a"] = "alta"
    mic_estimado_mg_l: float | None = None
    evidencias: list[ItemEvidenciaAPI] = Field(default_factory=list)
    mensagem: str | None = None


class RespostaPredicaoCompleta(BaseModel):
    """Contrato completo de resposta do job de predição (Seção 9.3, REG-14, REG-15)."""

    job_id: str
    status: Literal["em_fila", "processando", "concluido", "nao_suportado", "rejeitado_qc", "erro"] = (
        "concluido"
    )
    especie: str
    formato_entrada: Literal["FASTA", "FASTQ"] = "FASTA"
    qc: ResumoQCResposta
    predicoes: list[ItemPredicaoFarmacoAPI] = Field(default_factory=list)
    versao_modelo: str
    versoes_bancos: dict[str, str]
    latencia_ms: float | None = None
    aviso: str = AVISO_USO_PESQUISA


class InfoModeloCardAPI(BaseModel):
    """Informações e model card resumido expostos em GET /v1/modelos."""

    versao_modelo: str
    especie: str
    farmaco: str
    auprc_validacao_grupo: float
    vme_validacao_grupo: float
    me_validacao_grupo: float
    supera_baseline_b1: bool
    versoes_bancos: dict[str, str]
    aviso: str = AVISO_USO_PESQUISA


class RespostaObservabilidadeAPI(BaseModel):
    """Métricas de fila, latência (p50/p95) e monitoramento de drift (Seção 9.4)."""

    total_jobs_processados: int
    latencia_p50_ms: float
    latencia_p95_ms: float
    meta_latencia_p95_ms: float = 120_000.0
    dentro_meta_sla: bool = True
    distribuicao_especies: dict[str, int]
    taxa_aprovacao_qc: float
    probabilidade_resistencia_media: float
    taxa_abstencao_recente: float
    alerta_drift: bool = False
    aviso: str = AVISO_USO_PESQUISA
