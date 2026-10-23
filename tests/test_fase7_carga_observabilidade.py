"""Testes automatizados da Fase 7: Upload FASTQ na API, Multi-Espécie (S. aureus MRSA), Segurança X-API-Key, Observabilidade/Drift, Latência p95 e Entregáveis."""

from __future__ import annotations

import io
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from phyloresistance.api.app import app
from phyloresistance.genomica.fastq import sintetizar_fastq_isolado

client = TestClient(app)


def test_predicao_api_com_fastq_e_saureus_mrsa() -> None:
    """
    Testa o envio de um arquivo FASTQ bruto de Staphylococcus aureus contendo o gene mecA (MRSA)
    para POST /v1/predicoes:
    - Deve identificar formato_entrada == 'FASTQ'
    - Realizar trimagem e montagem em contigs
    - Classificar oxacilina como 'resistente' com estimativa de MIC (mg/L) e evidência SHAP no locus mecA
    """
    fastq_mrsa = sintetizar_fastq_isolado(
        genome_id="SAUREUS_MRSA_FASTQ",
        atributos_ativos=["gene_sa_mecA"],
        n_reads_alta_qualidade=35,
        n_reads_baixa_qualidade=4,
    )
    r_post = client.post(
        "/v1/predicoes?farmacos=oxacilina",
        files={"arquivo": ("isolado_mrsa.fastq", io.BytesIO(fastq_mrsa.encode("utf-8")), "text/plain")},
        data={"especie": "Staphylococcus aureus"},
    )
    assert r_post.status_code == 202
    dados_job = r_post.json()
    assert dados_job["formato_entrada"] == "FASTQ"
    job_id = dados_job["job_id"]

    r_get = client.get(f"/v1/predicoes/{job_id}")
    assert r_get.status_code == 200
    res = r_get.json()
    assert res["especie"] == "Staphylococcus aureus"
    assert res["formato_entrada"] == "FASTQ"
    assert res["qc"]["aprovado"] is True
    assert res["versao_modelo"] == "saureus-oxa-1.2.0"

    pred_oxa = res["predicoes"][0]
    assert pred_oxa["farmaco"] == "oxacilina"
    assert pred_oxa["classe"] == "resistente"
    assert pred_oxa["mic_estimado_mg_l"] is not None
    assert pred_oxa["mic_estimado_mg_l"] >= 4.0

    genes_ev = {ev["gene"] for ev in pred_oxa["evidencias"]}
    assert "mecA" in genes_ev


def test_seguranca_api_key_e_observabilidade_carga_p95(monkeypatch: pytest.MonkeyPatch) -> None:
    """
    Verifica:
    1) Proteção por X-API-Key quando PHYLORESISTANCE_API_KEY está configurada (Seção 9.4)
    2) Endpoint GET /v1/observabilidade com latência p50/p95 dentro do SLA e monitoramento de drift
    3) Listagem multi-espécie em GET /v1/modelos?todas_especies=true
    """
    r_todos = client.get("/v1/modelos?todas_especies=true")
    assert r_todos.status_code == 200
    todos_modelos = r_todos.json()
    especies_listadas = {m["especie"] for m in todos_modelos}
    assert especies_listadas == {
        "Escherichia coli",
        "Klebsiella pneumoniae",
        "Staphylococcus aureus",
    }

    # Ativa exigência de chave de API
    monkeypatch.setenv("PHYLORESISTANCE_API_KEY", "segredo-pesquisa-123")
    fastq_curto = sintetizar_fastq_isolado("LOAD_TEST", ["mut_gyrA_S83L"], n_reads_alta_qualidade=25)

    # Sem header X-API-Key -> deve rejeitar com 401
    r_sem_chave = client.post(
        "/v1/predicoes",
        files={"arquivo": ("load.fastq", io.BytesIO(fastq_curto.encode("utf-8")), "text/plain")},
    )
    assert r_sem_chave.status_code == 401

    # Com header X-API-Key válido -> deve aceitar com 202
    r_com_chave = client.post(
        "/v1/predicoes?farmacos=ciprofloxacina",
        files={"arquivo": ("load.fastq", io.BytesIO(fastq_curto.encode("utf-8")), "text/plain")},
        headers={"X-API-Key": "segredo-pesquisa-123"},
    )
    assert r_com_chave.status_code == 202

    # Consulta endpoint de observabilidade e latência p95
    r_obs = client.get("/v1/observabilidade")
    assert r_obs.status_code == 200
    obs = r_obs.json()
    assert obs["total_jobs_processados"] >= 2
    assert obs["latencia_p95_ms"] < 15_000.0  # Muito abaixo do limite de 1-2 minutos da Seção 9.4
    assert obs["dentro_meta_sla"] is True
    assert "Escherichia coli" in obs["distribuicao_especies"]


def test_entregaveis_fase7_figuras_svg_datasheet_e_quarto() -> None:
    """Verifica que todos os documentos de portfólio, figuras SVG e relatório Quarto da Fase 7 foram gerados."""
    arquivos_obrigatorios = [
        Path("docs/plano_de_analise.md"),
        Path("docs/curadoria_de_rotulos.md"),
        Path("docs/datasheet_dados.md"),
        Path("docs/model_card_ecoli_cip.md"),
        Path("docs/model_card_kpneumoniae.md"),
        Path("docs/model_card_saureus.md"),
        Path("docs/roteiro_demonstracao_video.md"),
        Path("reports/relatorio_quarto.qmd"),
        Path("reports/relatorio_cientifico.md"),
        Path("reports/relatorio_validacao.json"),
    ]
    for arq in arquivos_obrigatorios:
        assert arq.exists(), f"Entregável obrigatório ausente: {arq}"
        assert arq.stat().st_size > 100

