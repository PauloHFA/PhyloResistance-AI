"""Contratos de dados Pydantic para ingestão, curadoria de rótulos e versionamento por hash."""

from __future__ import annotations

from datetime import UTC, datetime
from typing import Literal

from pydantic import BaseModel, Field, field_validator


class RegistroBrutoBVBRC(BaseModel):
    """Registro bruto de genoma e teste de sensibilidade antimicrobiana (AST) vindo do BV-BRC/NCBI."""

    genome_id: str = Field(..., min_length=1, description="ID único do isolado bacteriano")
    especie: str = Field(..., min_length=1, description="Nome científico da espécie declarada")
    farmaco: str = Field(..., min_length=1, description="Nome do antimicrobiano testado")
    fenotipo_bruto: str | None = Field(
        default=None, description="Rótulo bruto (Resistant, Susceptible, Intermediate, R, S, I)"
    )
    mic_mg_l: float | None = Field(
        default=None, description="Concentração Inibitória Mínima (MIC) em mg/L"
    )
    metodo_laboratorio: str = Field(
        default="Broth microdilution",
        description="Método laboratorial (Broth microdilution, Disk diffusion, Agar dilution, Etest)",
    )
    padrao: str = Field(default="CLSI", description="Padrão de interpretação (CLSI ou EUCAST)")
    versao_ponto_corte: str = Field(
        default="2024.1", description="Versão do documento de ponto de corte"
    )
    pais: str = Field(default="Desconhecido", description="País de origem do isolado")
    ano: int = Field(default=2022, ge=1900, le=2100, description="Ano de coleta do isolado")
    st: str = Field(default="ST_ND", description="Sequence Type (MLST) do isolado")
    caminho_fasta: str | None = Field(
        default=None, description="Caminho opcional para o arquivo FASTA montado"
    )
    sha256_fasta: str | None = Field(
        default=None, description="Hash SHA-256 do arquivo FASTA bruto"
    )

    @field_validator("farmaco")
    @classmethod
    def normalizar_nome_farmaco(cls, v: str) -> str:
        mapa = {
            "ciprofloxacin": "ciprofloxacina",
            "ciprofloxacina": "ciprofloxacina",
            "cip": "ciprofloxacina",
            "gentamicin": "gentamicina",
            "gentamicina": "gentamicina",
            "gen": "gentamicina",
            "cefuroxime": "cefuroxima",
            "cefuroxima": "cefuroxima",
            "cxm": "cefuroxima",
            "ceftriaxone": "ceftriaxona",
            "ceftriaxona": "ceftriaxona",
            "cro": "ceftriaxona",
            "oxacillin": "oxacilina",
            "oxacilina": "oxacilina",
            "oxa": "oxacilina",
        }
        chave = v.strip().lower()
        return mapa.get(chave, chave)


class RegistroCurado(BaseModel):
    """Registro fenotípico padronizado e aprovado pela curadoria (REG-06)."""

    genome_id: str
    especie: str
    farmaco: str
    fenotipo: Literal["R", "S"]
    y: Literal[0, 1]
    origem_rotulo: Literal["recalculado_mic", "fenotipo_declarado"]
    mic_mg_l: float | None = None
    metodo_laboratorio: str
    padrao: str
    versao_ponto_corte: str
    pais: str
    ano: int
    st: str


class FunilExclusao(BaseModel):
    """Registro auditável do funil de exclusão de rótulos e isolados (Seção 3.3 e REG-06)."""

    total_registros_brutos: int = 0
    removidos_especie_fora_escopo: int = 0
    removidos_farmaco_fora_escopo: int = 0
    removidos_fenotipo_invalido: int = 0
    removidos_intermediario: int = 0
    removidos_duplicata_concordante: int = 0
    removidos_conflito_rotulo: int = 0
    removidos_qc_genoma: int = 0
    removidos_desduplicacao_clonal: int = 0
    total_aprovados_curadoria: int = 0
    total_aprovados_pos_qc_dedup: int = 0
    modo_intermediario: str = "excluir"
    detalhes_por_farmaco: dict[str, dict[str, int]] = Field(default_factory=dict)


class ManifestoVersaoDados(BaseModel):
    """Manifesto de proveniência, versão de bancos e hashes SHA-256 (Fase 0 e REG-13)."""

    gerado_em_utc: str = Field(
        default_factory=lambda: datetime.now(UTC).isoformat()
    )
    fonte_primaria: str = "BV-BRC / CARD / NCBI RefSeq"
    versoes_bancos: dict[str, str]
    arquivos_hashes_sha256: dict[str, str]
    termos_uso_confirmados: bool = True
    aviso_regulatorio: str = "Uso exclusivo em pesquisa científica. Não substitui o antibiograma."
