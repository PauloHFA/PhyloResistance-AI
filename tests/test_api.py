"""Testes dos endpoints da API FastAPI assíncrona (Seção 9, REG-14 e REG-15)."""

from __future__ import annotations

import io

from fastapi.testclient import TestClient

from phyloresistance.api.app import app
from phyloresistance.config import AVISO_USO_PESQUISA, VERSOES_BANCOS
from phyloresistance.genomica.anotacao import sintetizar_fasta_isolado

client = TestClient(app)


def test_health_e_ready_endpoints() -> None:
    """Verifica endpoints de vivacidade (/health) e prontidão (/ready)."""
    r_health = client.get("/health")
    assert r_health.status_code == 200
    assert r_health.json()["status"] == "ok"

    r_ready = client.get("/ready")
    assert r_ready.status_code == 200
    dados = r_ready.json()
    assert dados["status"] == "ready"
    assert "ciprofloxacina" in dados["modelos_carregados"]
    assert dados["versoes_bancos"] == VERSOES_BANCOS


def test_listar_modelos_e_model_cards() -> None:
    """Verifica GET /v1/modelos com versões e métricas de validação por grupo."""
    resp = client.get("/v1/modelos")
    assert resp.status_code == 200
    modelos = resp.json()
    assert len(modelos) == 3
    farmacos = {m["farmaco"] for m in modelos}
    assert farmacos == {"ciprofloxacina", "gentamicina", "cefuroxima"}
    for m in modelos:
        assert m["aviso"] == AVISO_USO_PESQUISA
        assert "amrfinderplus" in m["versoes_bancos"]


def test_fluxo_completo_predicao_fasta_e_relatorio_reg14_reg15() -> None:
    """
    Testa o fluxo assíncrono ponta a ponta (Seção 9.2 e 9.3, REG-14, REG-15):
    1) POST /v1/predicoes com FASTA contendo mutação gyrA S83L + parC S80I -> HTTP 202 + job_id
    2) GET /v1/predicoes/{job_id} -> resultado com probabilidade, evidência SHAP no locus gyrA e versões
    3) GET /v1/predicoes/{job_id}/relatorio -> relatório HTML
    """
    fasta_resistente = sintetizar_fasta_isolado(
        genome_id="ECOLI_API_CIP_R",
        atributos_ativos=["mut_gyrA_S83L", "mut_gyrA_D87N", "mut_parC_S80I", "gene_blaCTX_M_15"],
        n_contigs=35,
        tamanho_total_pb=4_950_000,
    )
    arquivo_bytes = io.BytesIO(fasta_resistente.encode("utf-8"))

    r_post = client.post(
        "/v1/predicoes?farmacos=ciprofloxacina",
        files={"arquivo": ("isolado_ecoli.fasta", arquivo_bytes, "text/plain")},
        data={"especie": "Escherichia coli"},
    )
    assert r_post.status_code == 202
    job_info = r_post.json()
    assert "job_id" in job_info
    job_id = job_info["job_id"]

    r_get = client.get(f"/v1/predicoes/{job_id}")
    assert r_get.status_code == 200
    resultado = r_get.json()

    assert resultado["job_id"] == job_id
    assert resultado["especie"] == "Escherichia coli"
    assert resultado["qc"]["aprovado"] is True
    assert resultado["versao_modelo"] == "ecoli-cip-1.2.0"
    assert resultado["versoes_bancos"] == VERSOES_BANCOS
    assert resultado["aviso"] == AVISO_USO_PESQUISA

    preds = resultado["predicoes"]
    assert len(preds) == 1
    pred_cip = preds[0]
    assert pred_cip["farmaco"] == "ciprofloxacina"
    assert pred_cip["classe"] == "resistente"
    assert pred_cip["probabilidade_resistencia"] >= 0.70

    genes_evidencia = {ev["gene"] for ev in pred_cip["evidencias"]}
    assert "gyrA" in genes_evidencia or "parC" in genes_evidencia

    # Verifica geração do relatório HTML
    r_html = client.get(f"/v1/predicoes/{job_id}/relatorio")
    assert r_html.status_code == 200
    assert "PhyloResistance AI" in r_html.text
    assert AVISO_USO_PESQUISA in r_html.text


def test_rejeicao_fasta_invalido_e_especie_fora_do_dominio() -> None:
    """Verifica rejeição HTTP 422 para arquivo não-FASTA e retorno 'nao_suportado' para espécie fora do domínio."""
    r_invalido = client.post(
        "/v1/predicoes",
        files={"arquivo": ("invalido.txt", io.BytesIO(b"isto nao eh um fasta"), "text/plain")},
    )
    assert r_invalido.status_code == 422

    fasta_ok = sintetizar_fasta_isolado(
        genome_id="ACINETO_01",
        atributos_ativos=[],
        n_contigs=30,
        tamanho_total_pb=4_800_000,
    )
    r_post = client.post(
        "/v1/predicoes",
        files={"arquivo": ("acineto.fasta", io.BytesIO(fasta_ok.encode("utf-8")), "text/plain")},
        data={"especie": "Acinetobacter baumannii"},
    )
    assert r_post.status_code == 202
    job_id = r_post.json()["job_id"]

    r_get = client.get(f"/v1/predicoes/{job_id}")
    assert r_get.status_code == 200
    dados = r_get.json()
    assert dados["status"] == "nao_suportado"
    assert dados["predicoes"][0]["classe"] == "nao_suportado"

