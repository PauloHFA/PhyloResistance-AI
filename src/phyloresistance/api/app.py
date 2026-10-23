"""Aplicação FastAPI assíncrona com fila de jobs, validação FASTA, QC e relatórios HTML/JSON (Seção 9)."""

from __future__ import annotations

import tempfile
import time
import uuid
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Annotated

from fastapi import FastAPI, File, Form, HTTPException, Query, UploadFile
from fastapi.responses import HTMLResponse

from phyloresistance.api.esquemas import (
    InfoModeloCardAPI,
    ItemEvidenciaAPI,
    ItemPredicaoFarmacoAPI,
    RespostaJobCriado,
    RespostaPredicaoCompleta,
    ResumoQCResposta,
)
from phyloresistance.config import (
    AVISO_USO_PESQUISA,
    ESPECIE_PILOTO,
    FARMACOS_PILOTO,
    VERSOES_BANCOS,
    VERSOES_MODELOS,
)
from phyloresistance.explicabilidade.shap_locus import gerar_relatorio_html_isolado
from phyloresistance.genomica.qc import (
    ErroValidacaoFASTA,
    avaliar_qc_genoma,
    validar_fasta_entrada,
)
from phyloresistance.pipeline import (
    RegistroModelosServico,
    executar_pipeline_completo,
)

LIMITE_MB_PADRAO: float = 25.0

_REGISTRO_GLOBAL: RegistroModelosServico | None = None
_ARMAZENAMENTO_JOBS: dict[str, dict[str, object]] = {}
_EXECUTOR_FILA = ThreadPoolExecutor(max_workers=4)


def obter_registro_modelos() -> RegistroModelosServico:
    """Carrega o registro de modelos de data/processed/modelos_servico.pkl ou treina coorte rápida sob demanda."""
    global _REGISTRO_GLOBAL
    if _REGISTRO_GLOBAL is not None and _REGISTRO_GLOBAL.esta_pronto():
        return _REGISTRO_GLOBAL

    caminho_padrao = Path("data/processed/modelos_servico.pkl")
    if caminho_padrao.exists():
        reg = RegistroModelosServico()
        reg.carregar(caminho_padrao)
        _REGISTRO_GLOBAL = reg
        return _REGISTRO_GLOBAL

    # Treinamento rápido de inicialização caso ainda não haja artefato em disco
    reg, _ = executar_pipeline_completo(
        output_dir="data/processed",
        reports_dir="reports",
        n_isolados_base=180,
        n_boot=40,
    )
    _REGISTRO_GLOBAL = reg
    return _REGISTRO_GLOBAL


async def salvar_temporario(arquivo: UploadFile, limite_mb: float = LIMITE_MB_PADRAO) -> Path:
    """Salva o upload em arquivo temporário isolado, verificando o limite máximo em MB."""
    conteudo = await arquivo.read()
    if not conteudo:
        raise ErroValidacaoFASTA("Arquivo enviado está vazio.")
    if len(conteudo) > limite_mb * 1024 * 1024:
        raise ErroValidacaoFASTA(f"Arquivo excede o limite de {limite_mb} MB.")

    tmp_dir = Path(tempfile.gettempdir()) / "phyloresistance_uploads"
    tmp_dir.mkdir(parents=True, exist_ok=True)
    caminho = tmp_dir / f"upload_{uuid.uuid4().hex[:12]}.fasta"
    caminho.write_bytes(conteudo)
    return caminho


def processar_genoma_worker(
    job_id: str,
    caminho_fasta: Path,
    especie: str,
    farmacos: list[str] | None,
) -> dict[str, object]:
    """
    Worker assíncrono que executa:
    1) Controle de Qualidade (QC) com os mesmos critérios do treino (Seção 9.4 e REG-07)
    2) Verificação de domínio de aplicação (espécie e fármacos suportados)
    3) Extração de atributos F1+F2+F3, predição calibrada com abstenção conformal e TreeSHAP agregado por locus
    4) Remoção imediata do arquivo temporário (Privacidade — Seção 9.4)
    """
    t0 = time.perf_counter()
    _ARMAZENAMENTO_JOBS[job_id] = {"job_id": job_id, "status": "processando"}

    try:
        stats = validar_fasta_entrada(caminho_fasta, limite_mb=LIMITE_MB_PADRAO)
        registro = obter_registro_modelos()

        # Verifica se a espécie está no domínio suportado
        especie_norm = especie.strip() if especie else ESPECIE_PILOTO
        qc_res = avaliar_qc_genoma(
            genome_id=job_id,
            tamanho_pb=int(stats["tamanho_pb"]),
            contigs=int(stats["contigs"]),
            n50=int(stats["n50"]),
            fracao_n=float(stats["fracao_n"]),
            especie_declarada=especie_norm,
            especie_inferida=especie_norm if especie_norm in registro.especies_suportadas else ESPECIE_PILOTO,
        )

        qc_payload = ResumoQCResposta(
            tamanho_pb=qc_res.tamanho_pb,
            contigs=qc_res.contigs,
            n50=qc_res.n50,
            aprovado=qc_res.aprovado,
            motivos_reprovacao=qc_res.motivos_reprovacao,
        )

        lista_farmacos = farmacos if farmacos else list(FARMACOS_PILOTO)
        versao_principal = (
            VERSOES_MODELOS.get(lista_farmacos[0], "ecoli-amr-1.2.0")
            if len(lista_farmacos) == 1
            else "ecoli-amr-1.2.0"
        )

        if especie_norm not in registro.especies_suportadas:
            resposta_fora_dominio = RespostaPredicaoCompleta(
                job_id=job_id,
                status="nao_suportado",
                especie=especie_norm,
                qc=qc_payload,
                predicoes=[
                    ItemPredicaoFarmacoAPI(
                        farmaco=f,
                        probabilidade_resistencia=None,
                        classe="nao_suportado",
                        confianca="n/a",
                        evidencias=[],
                        mensagem=f"Espécie '{especie_norm}' não suportada. Domínio válido: {sorted(registro.especies_suportadas)}.",
                    )
                    for f in lista_farmacos
                ],
                versao_modelo=versao_principal,
                versoes_bancos=dict(VERSOES_BANCOS),
                latencia_ms=round((time.perf_counter() - t0) * 1000.0, 2),
                aviso=AVISO_USO_PESQUISA,
            )
            payload = resposta_fora_dominio.model_dump()
            _ARMAZENAMENTO_JOBS[job_id] = payload
            return payload

        preds_brutas = registro.predizer_isolado_fasta(
            caminho_fasta, farmacos_solicitados=lista_farmacos
        )
        itens_pred: list[ItemPredicaoFarmacoAPI] = []
        for pb in preds_brutas:
            if pb.get("status") == "nao_suportado":
                itens_pred.append(
                    ItemPredicaoFarmacoAPI(
                        farmaco=str(pb["farmaco"]),
                        probabilidade_resistencia=None,
                        classe="nao_suportado",
                        confianca="n/a",
                        evidencias=[],
                        mensagem=str(pb.get("mensagem", "Fármaco não suportado.")),
                    )
                )
            else:
                evs_raw = pb.get("evidencias", [])
                evs_api = [
                    ItemEvidenciaAPI(**ev)  # type: ignore[arg-type]
                    for ev in (evs_raw if isinstance(evs_raw, list) else [])
                ]
                itens_pred.append(
                    ItemPredicaoFarmacoAPI(
                        farmaco=str(pb["farmaco"]),
                        probabilidade_resistencia=float(pb["probabilidade_resistencia"]),  # type: ignore[arg-type]
                        classe=str(pb["classe"]),  # type: ignore[arg-type]
                        confianca=str(pb["confianca"]),  # type: ignore[arg-type]
                        evidencias=evs_api,
                    )
                )

        status_final = "concluido" if qc_res.aprovado else "rejeitado_qc"
        resposta = RespostaPredicaoCompleta(
            job_id=job_id,
            status=status_final,  # type: ignore[arg-type]
            especie=especie_norm,
            qc=qc_payload,
            predicoes=itens_pred,
            versao_modelo=versao_principal,
            versoes_bancos=dict(VERSOES_BANCOS),
            latencia_ms=round((time.perf_counter() - t0) * 1000.0, 2),
            aviso=AVISO_USO_PESQUISA,
        )
        payload = resposta.model_dump()
        _ARMAZENAMENTO_JOBS[job_id] = payload
        return payload
    finally:
        # Privacidade (Seção 9.4): apaga imediatamente o arquivo FASTA temporário
        try:
            if caminho_fasta.exists():
                caminho_fasta.unlink()
        except OSError:
            pass


def criar_app() -> FastAPI:
    """Cria e configura a aplicação FastAPI do PhyloResistance AI."""
    aplicacao = FastAPI(
        title="PhyloResistance AI API",
        version="1.2.0",
        description=(
            "API assíncrona para predição de Resistência Antimicrobiana (AMR) e "
            "mapeamento de evidências genômicas por Machine Learning Explicável (TreeSHAP). "
            f"{AVISO_USO_PESQUISA}"
        ),
    )

    @aplicacao.get("/health")
    async def health() -> dict[str, str]:
        """Vivacidade do serviço (Seção 9.2)."""
        return {"status": "ok", "servico": "phyloresistance-api"}

    @aplicacao.get("/ready")
    async def ready() -> dict[str, object]:
        """Prontidão do serviço: verifica se modelos e catálogos estão carregados (Seção 9.2)."""
        reg = obter_registro_modelos()
        return {
            "status": "ready" if reg.esta_pronto() else "not_ready",
            "modelos_carregados": list(reg.modelos_por_farmaco.keys()),
            "versoes_bancos": dict(VERSOES_BANCOS),
            "aviso": AVISO_USO_PESQUISA,
        }

    @aplicacao.get("/v1/modelos", response_model=list[InfoModeloCardAPI])
    async def listar_modelos() -> list[InfoModeloCardAPI]:
        """Lista modelos disponíveis, versões, métricas de validação por grupo e model cards (Seção 9.2)."""
        reg = obter_registro_modelos()
        saida: list[InfoModeloCardAPI] = []
        for farmaco, art in reg.modelos_por_farmaco.items():
            saida.append(
                InfoModeloCardAPI(
                    versao_modelo=art.versao_modelo,
                    especie=art.especie,
                    farmaco=farmaco,
                    auprc_validacao_grupo=art.metricas_validacao_principal.auprc,
                    vme_validacao_grupo=art.metricas_validacao_principal.vme,
                    me_validacao_grupo=art.metricas_validacao_principal.me,
                    supera_baseline_b1=art.comparacao_b1.supera_b1_com_ic95,
                    versoes_bancos=dict(VERSOES_BANCOS),
                    aviso=AVISO_USO_PESQUISA,
                )
            )
        return saida

    @aplicacao.post("/v1/predicoes", status_code=202, response_model=RespostaJobCriado)
    async def criar_predicao(
        arquivo: Annotated[UploadFile, File(description="Genoma bacteriano montado em formato FASTA")],
        especie: Annotated[str | None, Form()] = None,
        farmacos: Annotated[list[str] | None, Query()] = None,
    ) -> RespostaJobCriado:
        """
        Envia um arquivo FASTA montado (e opcionalmente espécie e lista de fármacos);
        valida o formato e alfabeto ACGTN e enfileira o job assíncrono retornando 202 com job_id (Seção 9.2).
        """
        caminho_tmp = await salvar_temporario(arquivo, limite_mb=LIMITE_MB_PADRAO)
        try:
            validar_fasta_entrada(caminho_tmp, limite_mb=LIMITE_MB_PADRAO)
        except ErroValidacaoFASTA as exc:
            if caminho_tmp.exists():
                caminho_tmp.unlink()
            raise HTTPException(status_code=422, detail=str(exc)) from exc

        job_id = uuid.uuid4().hex[:8]
        _ARMAZENAMENTO_JOBS[job_id] = {"job_id": job_id, "status": "em_fila"}

        # Submete à fila de workers
        future = _EXECUTOR_FILA.submit(
            processar_genoma_worker,
            job_id,
            caminho_tmp,
            especie or ESPECIE_PILOTO,
            farmacos,
        )
        # Aguarda brevemente caso o worker local já conclua em poucos milissegundos
        try:
            future.result(timeout=2.0)
        except Exception:
            pass

        return RespostaJobCriado(job_id=job_id, status="em_fila")

    @aplicacao.get("/v1/predicoes/{job_id}")
    async def consultar_predicao(job_id: str) -> dict[str, object]:
        """Consulta o estado do job e, quando pronto, retorna o resultado completo (Seção 9.2 e 9.3)."""
        if job_id not in _ARMAZENAMENTO_JOBS:
            raise HTTPException(status_code=404, detail=f"Job '{job_id}' não encontrado.")
        return _ARMAZENAMENTO_JOBS[job_id]

    @aplicacao.get("/v1/predicoes/{job_id}/relatorio", response_class=HTMLResponse)
    async def consultar_relatorio_predicao(job_id: str) -> HTMLResponse:
        """Gera e retorna o relatório interpretável em HTML para o isolado (Seção 9.2)."""
        if job_id not in _ARMAZENAMENTO_JOBS:
            raise HTTPException(status_code=404, detail=f"Job '{job_id}' não encontrado.")
        dados = _ARMAZENAMENTO_JOBS[job_id]
        html = gerar_relatorio_html_isolado(dados)
        return HTMLResponse(content=html, status_code=200)

    return aplicacao


app = criar_app()

