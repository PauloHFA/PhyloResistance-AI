"""Aplicação FastAPI assíncrona com fila de jobs, suporte FASTA/FASTQ, QC, observabilidade/drift e segurança (Seção 9)."""

from __future__ import annotations

import os
import tempfile
import time
import uuid
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Annotated, Literal

import numpy as np
from fastapi import FastAPI, File, Form, Header, HTTPException, Query, Request, UploadFile
from fastapi.responses import HTMLResponse, Response

from phyloresistance.api.esquemas import (
    InfoModeloCardAPI,
    ItemEvidenciaAPI,
    ItemPredicaoFarmacoAPI,
    RespostaJobCriado,
    RespostaObservabilidadeAPI,
    RespostaPredicaoCompleta,
    ResumoQCResposta,
)
from phyloresistance.config import (
    AVISO_USO_PESQUISA,
    CRITERIOS_QC_POR_ESPECIE,
    ESPECIE_PILOTO,
    FARMACOS_POR_ESPECIE,
    VERSOES_BANCOS,
    VERSOES_MODELOS_POR_ESPECIE,
)
from phyloresistance.explicabilidade.shap_locus import (
    gerar_relatorio_html_isolado,
    gerar_relatorio_pdf_isolado,
)
from phyloresistance.genomica.fastq import (
    ErroValidacaoFASTQ,
    trimar_e_montar_fastq_para_fasta,
)
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
LIMITE_RATE_POR_MINUTO: int = 120

_REGISTRO_GLOBAL: RegistroModelosServico | None = None
_ARMAZENAMENTO_JOBS: dict[str, dict[str, object]] = {}
_EXECUTOR_FILA = ThreadPoolExecutor(max_workers=4)

# Estado de observabilidade e monitoramento de drift (Seção 9.4)
_LATENCIAS_MS: list[float] = []
_CONTAGEM_ESPECIES: dict[str, int] = defaultdict(int)
_HISTORICO_QC_APROVADO: list[int] = []
_HISTORICO_PROBS: list[float] = []
_HISTORICO_ABSTENCAO: list[int] = []
_JANELA_RATE_LIMIT: dict[str, list[float]] = defaultdict(list)


def verificar_seguranca_e_rate_limit(
    request: Request,
    x_api_key: str | None = None,
) -> None:
    """
    Verifica chave de API (quando PHYLORESISTANCE_API_KEY estiver configurada no ambiente)
    e aplica rate limit por IP/cliente (Seção 9.4).
    """
    chave_esperada = os.environ.get("PHYLORESISTANCE_API_KEY")
    if chave_esperada and x_api_key != chave_esperada:
        raise HTTPException(status_code=401, detail="Chave de API (X-API-Key) inválida ou ausente.")

    cliente_id = request.client.host if request.client else "local"
    agora = time.time()
    timestamps = [t for t in _JANELA_RATE_LIMIT[cliente_id] if (agora - t) < 60.0]
    if len(timestamps) >= LIMITE_RATE_POR_MINUTO:
        raise HTTPException(status_code=429, detail="Limite de requisições por minuto excedido (rate limit).")
    timestamps.append(agora)
    _JANELA_RATE_LIMIT[cliente_id] = timestamps


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

    reg, _ = executar_pipeline_completo(
        output_dir="data/processed",
        reports_dir="reports",
        n_isolados_base=150,
        n_boot=40,
    )
    _REGISTRO_GLOBAL = reg
    return _REGISTRO_GLOBAL


async def salvar_temporario(
    arquivo: UploadFile, limite_mb: float = LIMITE_MB_PADRAO
) -> tuple[Path, Literal["FASTA", "FASTQ"]]:
    """Salva o upload em arquivo temporário isolado e identifica se é FASTA ou FASTQ."""
    conteudo = await arquivo.read()
    if not conteudo:
        raise ErroValidacaoFASTA("Arquivo enviado está vazio.")
    if len(conteudo) > limite_mb * 1024 * 1024:
        raise ErroValidacaoFASTA(f"Arquivo excede o limite de {limite_mb} MB.")

    texto_ini = conteudo[:256].decode("utf-8", errors="ignore").lstrip()
    nome_arq = (arquivo.filename or "").lower()
    eh_fastq = texto_ini.startswith("@") or nome_arq.endswith(".fastq") or nome_arq.endswith(".fq")
    formato: Literal["FASTA", "FASTQ"] = "FASTQ" if eh_fastq else "FASTA"

    tmp_dir = Path(tempfile.gettempdir()) / "phyloresistance_uploads"
    tmp_dir.mkdir(parents=True, exist_ok=True)
    ext = ".fastq" if eh_fastq else ".fasta"
    caminho = tmp_dir / f"upload_{uuid.uuid4().hex[:12]}{ext}"
    caminho.write_bytes(conteudo)
    return caminho, formato


def processar_genoma_worker(
    job_id: str,
    caminho_entrada: Path,
    formato_entrada: Literal["FASTA", "FASTQ"],
    especie: str,
    farmacos: list[str] | None,
) -> dict[str, object]:
    """
    Worker assíncrono que executa:
    0) Se FASTQ, executa trimagem Phred e montagem em contigs FASTA (Extensão Fase 6)
    1) Controle de Qualidade (QC) com os mesmos critérios do treino (Seção 9.4 e REG-07)
    2) Verificação de domínio de aplicação (espécie e fármacos suportados)
    3) Extração de atributos F1+F2+F3, predição calibrada com abstenção conformal, MIC e TreeSHAP
    4) Registro de métricas de observabilidade/drift e remoção imediata do arquivo temporário
    """
    t0 = time.perf_counter()
    _ARMAZENAMENTO_JOBS[job_id] = {"job_id": job_id, "status": "processando"}
    caminho_fasta_gerado: Path | None = None

    try:
        especie_norm = especie.strip() if especie else ESPECIE_PILOTO
        crit_esp = CRITERIOS_QC_POR_ESPECIE.get(especie_norm)
        tam_alvo = (
            (crit_esp.tamanho_min_pb + crit_esp.tamanho_max_pb) // 2
            if crit_esp is not None
            else 4_850_000
        )

        if formato_entrada == "FASTQ":
            res_fq = trimar_e_montar_fastq_para_fasta(
                caminho_entrada,
                genome_id=job_id,
                tamanho_alvo_montagem_pb=tam_alvo,
            )
            caminho_fasta_gerado = caminho_entrada.with_suffix(".assembled.fasta")
            caminho_fasta_gerado.write_text(res_fq.fasta_montado, encoding="utf-8")
            caminho_para_qc = caminho_fasta_gerado
        else:
            caminho_para_qc = caminho_entrada

        stats = validar_fasta_entrada(caminho_para_qc, limite_mb=LIMITE_MB_PADRAO)
        registro = obter_registro_modelos()

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

        farmacos_padrao_esp = list(FARMACOS_POR_ESPECIE.get(especie_norm, ("ciprofloxacina",)))
        lista_farmacos = farmacos if farmacos else farmacos_padrao_esp
        mapa_versoes_esp = VERSOES_MODELOS_POR_ESPECIE.get(especie_norm, {})
        versao_principal = (
            mapa_versoes_esp.get(lista_farmacos[0], "ecoli-amr-1.2.0")
            if len(lista_farmacos) == 1
            else "phyloresistance-suite-1.2.0"
        )

        _CONTAGEM_ESPECIES[especie_norm] += 1
        _HISTORICO_QC_APROVADO.append(1 if qc_res.aprovado else 0)

        if especie_norm not in registro.especies_suportadas:
            lat = round((time.perf_counter() - t0) * 1000.0, 2)
            _LATENCIAS_MS.append(lat)
            resposta_fora_dominio = RespostaPredicaoCompleta(
                job_id=job_id,
                status="nao_suportado",
                especie=especie_norm,
                formato_entrada=formato_entrada,
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
                latencia_ms=lat,
                aviso=AVISO_USO_PESQUISA,
            )
            payload = resposta_fora_dominio.model_dump()
            _ARMAZENAMENTO_JOBS[job_id] = payload
            return payload

        preds_brutas = registro.predizer_isolado_fasta(
            caminho_para_qc,
            farmacos_solicitados=lista_farmacos,
            especie=especie_norm,
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
                prob_val = float(pb["probabilidade_resistencia"])  # type: ignore[arg-type]
                classe_val = str(pb["classe"])
                mic_est = pb.get("mic_estimado_mg_l")
                _HISTORICO_PROBS.append(prob_val)
                _HISTORICO_ABSTENCAO.append(1 if classe_val == "indeterminado" else 0)

                itens_pred.append(
                    ItemPredicaoFarmacoAPI(
                        farmaco=str(pb["farmaco"]),
                        probabilidade_resistencia=prob_val,
                        classe=classe_val,  # type: ignore[arg-type]
                        confianca=str(pb["confianca"]),  # type: ignore[arg-type]
                        mic_estimado_mg_l=float(mic_est) if isinstance(mic_est, (int, float)) else None,
                        evidencias=evs_api,
                    )
                )

        lat = round((time.perf_counter() - t0) * 1000.0, 2)
        _LATENCIAS_MS.append(lat)
        status_final = "concluido" if qc_res.aprovado else "rejeitado_qc"
        resposta = RespostaPredicaoCompleta(
            job_id=job_id,
            status=status_final,  # type: ignore[arg-type]
            especie=especie_norm,
            formato_entrada=formato_entrada,
            qc=qc_payload,
            predicoes=itens_pred,
            versao_modelo=versao_principal,
            versoes_bancos=dict(VERSOES_BANCOS),
            latencia_ms=lat,
            aviso=AVISO_USO_PESQUISA,
        )
        payload = resposta.model_dump()
        _ARMAZENAMENTO_JOBS[job_id] = payload
        return payload
    finally:
        # Privacidade (Seção 9.4): apaga imediatamente os arquivos temporários
        for cam in (caminho_entrada, caminho_fasta_gerado):
            try:
                if cam is not None and cam.exists():
                    cam.unlink()
            except OSError:
                pass


def criar_app() -> FastAPI:
    """Cria e configura a aplicação FastAPI do PhyloResistance AI."""
    aplicacao = FastAPI(
        title="PhyloResistance AI API",
        version="1.2.0",
        description=(
            "API assíncrona para predição de Resistência Antimicrobiana (AMR), estimativa de MIC e "
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
            "especies_suportadas": sorted(reg.especies_suportadas),
            "modelos_carregados": list(reg.modelos_por_farmaco.keys()),
            "total_modelos_multiespecie": len(reg.listar_todos_artefatos()),
            "versoes_bancos": dict(VERSOES_BANCOS),
            "aviso": AVISO_USO_PESQUISA,
        }

    @aplicacao.get("/v1/modelos", response_model=list[InfoModeloCardAPI])
    async def listar_modelos(
        especie: str | None = None,
        todas_especies: bool = False,
    ) -> list[InfoModeloCardAPI]:
        """Lista modelos disponíveis, versões, métricas de validação por grupo e model cards (Seção 9.2)."""
        reg = obter_registro_modelos()
        if todas_especies:
            artefatos = reg.listar_todos_artefatos()
        elif especie is not None:
            artefatos = [a for a in reg.listar_todos_artefatos() if a.especie == especie]
        else:
            artefatos = list(reg.modelos_por_farmaco.values())

        saida: list[InfoModeloCardAPI] = []
        for art in artefatos:
            saida.append(
                InfoModeloCardAPI(
                    versao_modelo=art.versao_modelo,
                    especie=art.especie,
                    farmaco=art.farmaco,
                    auprc_validacao_grupo=art.metricas_validacao_principal.auprc,
                    vme_validacao_grupo=art.metricas_validacao_principal.vme,
                    me_validacao_grupo=art.metricas_validacao_principal.me,
                    supera_baseline_b1=art.comparacao_b1.supera_b1_com_ic95,
                    versoes_bancos=dict(VERSOES_BANCOS),
                    aviso=AVISO_USO_PESQUISA,
                )
            )
        return saida

    @aplicacao.get("/v1/observabilidade", response_model=RespostaObservabilidadeAPI)
    async def consultar_observabilidade() -> RespostaObservabilidadeAPI:
        """
        Exibe métricas de latência (p50, p95), taxa de aprovação de QC, distribuição de espécies
        e monitoramento de drift de probabilidades/abstenção (Seção 9.4 e Fase 5/7).
        """
        n_jobs = len(_LATENCIAS_MS)
        p50 = float(np.percentile(_LATENCIAS_MS, 50)) if n_jobs > 0 else 0.0
        p95 = float(np.percentile(_LATENCIAS_MS, 95)) if n_jobs > 0 else 0.0
        taxa_qc = float(np.mean(_HISTORICO_QC_APROVADO)) if _HISTORICO_QC_APROVADO else 1.0
        prob_m = float(np.mean(_HISTORICO_PROBS)) if _HISTORICO_PROBS else 0.5
        taxa_abst = float(np.mean(_HISTORICO_ABSTENCAO)) if _HISTORICO_ABSTENCAO else 0.0

        # Sinaliza alerta de drift se a taxa de reprovação de QC > 50% ou taxa de abstenção > 45%
        alerta_drift = (taxa_qc < 0.50) or (taxa_abst > 0.45)

        return RespostaObservabilidadeAPI(
            total_jobs_processados=n_jobs,
            latencia_p50_ms=round(p50, 2),
            latencia_p95_ms=round(p95, 2),
            dentro_meta_sla=(p95 <= 120_000.0),
            distribuicao_especies=dict(_CONTAGEM_ESPECIES),
            taxa_aprovacao_qc=round(taxa_qc, 4),
            probabilidade_resistencia_media=round(prob_m, 4),
            taxa_abstencao_recente=round(taxa_abst, 4),
            alerta_drift=alerta_drift,
            aviso=AVISO_USO_PESQUISA,
        )

    @aplicacao.post("/v1/predicoes", status_code=202, response_model=RespostaJobCriado)
    async def criar_predicao(
        request: Request,
        arquivo: Annotated[UploadFile, File(description="Genoma bacteriano em formato FASTA ou FASTQ")],
        especie: Annotated[str | None, Form()] = None,
        farmacos: Annotated[list[str] | None, Query()] = None,
        x_api_key: Annotated[str | None, Header()] = None,
    ) -> RespostaJobCriado:
        """
        Envia um arquivo FASTA montado ou FASTQ bruto (e opcionalmente espécie e lista de fármacos);
        valida o formato e enfileira o job assíncrono retornando 202 com job_id (Seção 9.2 e Fase 6).
        """
        verificar_seguranca_e_rate_limit(request, x_api_key=x_api_key)
        caminho_tmp, formato = await salvar_temporario(arquivo, limite_mb=LIMITE_MB_PADRAO)
        try:
            if formato == "FASTQ":
                from phyloresistance.genomica.fastq import ler_registros_fastq

                ler_registros_fastq(caminho_tmp)
            else:
                validar_fasta_entrada(caminho_tmp, limite_mb=LIMITE_MB_PADRAO)
        except (ErroValidacaoFASTA, ErroValidacaoFASTQ) as exc:
            if caminho_tmp.exists():
                caminho_tmp.unlink()
            raise HTTPException(status_code=422, detail=str(exc)) from exc

        job_id = uuid.uuid4().hex[:8]
        _ARMAZENAMENTO_JOBS[job_id] = {
            "job_id": job_id,
            "status": "em_fila",
            "formato_entrada": formato,
        }

        future = _EXECUTOR_FILA.submit(
            processar_genoma_worker,
            job_id,
            caminho_tmp,
            formato,
            especie or ESPECIE_PILOTO,
            farmacos,
        )
        try:
            future.result(timeout=3.0)
        except Exception:
            pass

        return RespostaJobCriado(job_id=job_id, status="em_fila", formato_entrada=formato)

    @aplicacao.get("/v1/predicoes/{job_id}")
    async def consultar_predicao(job_id: str) -> dict[str, object]:
        """Consulta o estado do job e, quando pronto, retorna o resultado completo (Seção 9.2 e 9.3)."""
        if job_id not in _ARMAZENAMENTO_JOBS:
            raise HTTPException(status_code=404, detail=f"Job '{job_id}' não encontrado.")
        return _ARMAZENAMENTO_JOBS[job_id]

    @aplicacao.get("/v1/predicoes/{job_id}/relatorio")
    async def consultar_relatorio_predicao(
        job_id: str,
        formato: Literal["html", "pdf"] = "html",
    ) -> Response:
        """Gera e retorna o relatório interpretável em HTML ou PDF para o isolado (Seção 9.2)."""
        if job_id not in _ARMAZENAMENTO_JOBS:
            raise HTTPException(status_code=404, detail=f"Job '{job_id}' não encontrado.")
        dados = _ARMAZENAMENTO_JOBS[job_id]
        if formato == "pdf":
            pdf_bytes = gerar_relatorio_pdf_isolado(dados)
            return Response(
                content=pdf_bytes,
                media_type="application/pdf",
                headers={"Content-Disposition": f'inline; filename="relatorio_{job_id}.pdf"'},
                status_code=200,
            )
        html = gerar_relatorio_html_isolado(dados)
        return HTMLResponse(content=html, status_code=200)

    return aplicacao


app = criar_app()
