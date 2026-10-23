"""
Testes automatizados dos complementos avançados da stack (Seções 3.1, 4.1, 5.1, 8.1, 9.2, 12 e 16):
- Motor alternativo LightGBM em B3 e B4
- Mapeamento físico de k-mers/unitigs contra coordenadas genômicas de referência GFF3 (Bakta/Prokka)
- Geração de relatório em PDF 1.4 na API (`GET /v1/predicoes/{job_id}/relatorio?formato=pdf`)
- Gerenciador de fila de workers (RQ / ThreadPool)
- Rastreador de experimentos MLflow (`RastreadorExperimentosMLflow`)
- Cliente de download BV-BRC (`baixar_dados_bvbrc_api`) e comando `phyloresistance demo`
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import scipy.sparse as sp
from fastapi.testclient import TestClient

from phyloresistance.api.app import app
from phyloresistance.api.worker import GerenciadorFilaWorkers
from phyloresistance.cli import main as cli_main
from phyloresistance.dados.bvbrc import (
    ASSINATURAS_NUCLEOTIDICAS_31MER,
    baixar_dados_bvbrc_api,
)
from phyloresistance.explicabilidade.mapeamento_genomico import (
    alinhar_kmers_contra_referencia_gff3,
    carregar_anotacao_gff3,
    gerar_referencia_e_gff3_sintetico_ecoli,
)
from phyloresistance.genomica.anotacao import sintetizar_fasta_isolado
from phyloresistance.modelos.escada import ModeloArvoreCuradoB3, PipelineCompletoB4
from phyloresistance.relatorios_visuais import gerar_todos_documentos_e_dashboards
from phyloresistance.validacao.experimentos import RastreadorExperimentosMLflow


def test_motor_lightgbm_em_b3_e_b4() -> None:
    """Valida que os modelos B3 e B4 suportam LightGBM além de XGBoost (Seção 5.1 e 12)."""
    rng = np.random.default_rng(42)
    X_f12 = rng.integers(0, 2, size=(60, 8)).astype(np.float64)
    y = ((X_f12[:, 0] + X_f12[:, 1]) > 0).astype(np.int32)
    X_f3 = sp.csr_matrix(rng.integers(0, 2, size=(60, 20), dtype=np.int8))

    b3_lgbm = ModeloArvoreCuradoB3(motor="lightgbm", n_estimators=25, seed=42).fit(X_f12, y)
    probs_b3 = b3_lgbm.predict_proba(X_f12)[:, 1]
    assert probs_b3.shape == (60,)
    assert float(np.mean((probs_b3 >= 0.5) == y)) >= 0.80

    nomes_f12 = [f"gene_{i}" for i in range(8)]
    nomes_f3 = [f"unitig_k31_{i:04d}" for i in range(20)]
    b4_lgbm = PipelineCompletoB4(
        nomes_f1_f2=nomes_f12,
        nomes_f3=nomes_f3,
        max_kmers_f3=10,
        motor="lightgbm",
        seed=42,
    ).fit(X_f12, X_f3, y)
    probs_b4 = b4_lgbm.predict_proba(X_f12, X_f3)[:, 1]
    assert probs_b4.shape == (60,)
    assert float(np.mean((probs_b4 >= 0.5) == y)) >= 0.80


def test_mapeamento_kmers_coordenadas_referencia_gff3(tmp_path: Path) -> None:
    """Valida o parser GFF3 (Bakta/Prokka) e o alinhamento físico de k-mers contra o genoma de referência (Seção 4.1 e 8.1)."""
    fasta_ref, gff3_str = gerar_referencia_e_gff3_sintetico_ecoli()
    caminho_gff = tmp_path / "referencia_ecoli.gff3"
    caminho_gff.write_text(gff3_str, encoding="utf-8")
    anotacoes = carregar_anotacao_gff3(caminho_gff)
    assert len(anotacoes) >= 6
    assert any(a.gene == "gyrA" for a in anotacoes)

    mapeamentos = alinhar_kmers_contra_referencia_gff3(
        kmers_por_id={
            "unitig_gyrA": ASSINATURAS_NUCLEOTIDICAS_31MER["mut_gyrA_S83L"],
            "unitig_marR": ASSINATURAS_NUCLEOTIDICAS_31MER["unitig_cand_marR_G103S"],
            "unitig_plasmidial": "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT",
        },
        fasta_referencia=fasta_ref,
        gff3_referencia=gff3_str,
    )
    assert len(mapeamentos) == 3
    assert mapeamentos[0].gene_anotado == "gyrA"
    assert mapeamentos[0].identidade == 100.0
    assert mapeamentos[1].gene_anotado == "marR"
    assert mapeamentos[2].seqid == "contig_acessorio_plasmidial"


def test_relatorio_pdf_api_e_gerenciador_worker(tmp_path: Path) -> None:
    """Valida a exportação do relatório em formato PDF 1.4 na API e o GerenciadorFilaWorkers (Seção 9.2)."""
    cliente = TestClient(app)
    fasta_path = tmp_path / "isolado_pdf.fasta"
    fasta_str = sintetizar_fasta_isolado(
        genome_id="ISO_PDF_01",
        atributos_ativos=["mut_gyrA_S83L", "mut_parC_S80I"],
    )
    fasta_path.write_text(fasta_str, encoding="utf-8")
    with fasta_path.open("rb") as f:
        r_post = cliente.post(
            "/v1/predicoes",
            files={"arquivo": ("isolado_pdf.fasta", f, "text/plain")},
        )
    assert r_post.status_code == 202
    job_id = r_post.json()["job_id"]

    r_pdf = cliente.get(f"/v1/predicoes/{job_id}/relatorio?formato=pdf")
    assert r_pdf.status_code == 200
    assert r_pdf.headers["content-type"] == "application/pdf"
    assert r_pdf.content.startswith(b"%PDF-1.4")
    assert b"PhyloResistance AI" in r_pdf.content
    assert b"%%EOF" in r_pdf.content

    gerenciador = GerenciadorFilaWorkers(redis_url=None)
    assert gerenciador.modo_execucao == "threadpool_local"


def test_rastreador_mlflow_cliente_bvbrc_e_cli_demo(tmp_path: Path) -> None:
    """Valida o rastreador de experimentos MLflow, o cliente BV-BRC e o subcomando `demo` da CLI."""
    rastreador = RastreadorExperimentosMLflow(diretorio_mlruns=tmp_path / "mlruns")
    run = rastreador.registrar_execucao_farmaco(
        especie="Escherichia coli",
        farmaco="ciprofloxacina",
        parametros={"modelo": "B4", "motor": "xgboost"},
        metricas={"auprc": 0.992, "vme": 0.015},
    )
    assert (tmp_path / "mlruns" / "indice_experimentos.json").exists()
    assert (tmp_path / "mlruns" / rastreador.nome_experimento / run.run_id / "meta.json").exists()

    regs_dl, cam_dl, hash_dl = baixar_dados_bvbrc_api(
        especie="Escherichia coli",
        limite=15,
        diretorio_raw=tmp_path / "raw_bvbrc",
        timeout_s=1.0,
    )
    assert len(regs_dl) == 15
    assert cam_dl.exists()
    assert len(hash_dl) == 64

    cod_saida = cli_main(["demo", "--especie", "Escherichia coli", "--output-dir", str(tmp_path)])
    assert cod_saida == 0


def test_geracao_documento_escrito_dashboard_html_e_pdf_multipagina() -> None:
    """Valida que o documento escrito (.md), o dashboard interativo (.html) e o relatório multipágina (.pdf) são gerados."""
    caminhos = gerar_todos_documentos_e_dashboards()
    cam_md = Path(caminhos["documento_escrito_md"])
    cam_html = Path(caminhos["dashboard_html_docs"])
    cam_pdf = Path(caminhos["relatorio_pdf_docs"])

    assert cam_md.exists() and cam_md.stat().st_size > 8_000
    assert cam_html.exists() and cam_html.stat().st_size > 20_000
    html_txt = cam_html.read_text(encoding="utf-8")
    assert "<svg" in html_txt and "Mapa Epidemiológico Global" in html_txt
    assert "NC_000913.3" in html_txt

    assert cam_pdf.exists() and cam_pdf.stat().st_size > 10_000
    pdf_bytes = cam_pdf.read_bytes()
    assert pdf_bytes.startswith(b"%PDF-1.4")
    assert b"/Count 4" in pdf_bytes
    assert b"%%EOF" in pdf_bytes

    cam_estudo_md = Path(caminhos["estudo_cientifico_md"])
    cam_estudo_html = Path(caminhos["estudo_cientifico_html"])
    cam_estudo_pdf = Path(caminhos["estudo_cientifico_pdf"])

    assert cam_estudo_md.exists() and cam_estudo_md.stat().st_size > 15_000
    txt_estudo_md = cam_estudo_md.read_text(encoding="utf-8")
    assert "Estudo Científico" in txt_estudo_md
    assert "Monografia" not in txt_estudo_md
    assert "Conclusão Geral: Onde Chegamos, Problemas Enfrentados" in txt_estudo_md
    assert "Antimicrobial Resistance Collaborators" in txt_estudo_md
    assert cam_estudo_html.exists() and cam_estudo_html.stat().st_size > 25_000
    assert cam_estudo_pdf.exists() and cam_estudo_pdf.stat().st_size > 12_000
    pdf_estudo_bytes = cam_estudo_pdf.read_bytes()
    assert pdf_estudo_bytes.startswith(b"%PDF-1.4")
    assert b"/Count 6" in pdf_estudo_bytes

