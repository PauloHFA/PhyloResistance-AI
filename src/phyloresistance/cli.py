"""Interface de Linha de Comando (CLI) para reprodução do pipeline PhyloResistance AI."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from phyloresistance.dados.bvbrc import baixar_dados_bvbrc_api
from phyloresistance.genomica.anotacao import sintetizar_fasta_isolado
from phyloresistance.genomica.qc import avaliar_qc_genoma, validar_fasta_entrada
from phyloresistance.pipeline import RegistroModelosServico, executar_pipeline_completo


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="phyloresistance",
        description="PhyloResistance AI — Pipeline reproduzível de predição de AMR e XAI genômico.",
    )
    sub = parser.add_subparsers(dest="comando")

    p_run = sub.add_parser("run-all", help="Executa o pipeline completo ponta a ponta (Fases 1 a 7).")
    p_run.add_argument("--output-dir", default="data/processed", help="Diretório de saída dos artefatos")
    p_run.add_argument("--reports-dir", default="reports", help="Diretório de saída dos relatórios")
    p_run.add_argument("--n-isolados", type=int, default=360, help="Tamanho da coorte piloto")
    p_run.add_argument("--n-boot", type=int, default=200, help="Iterações de bootstrap em blocos")

    p_demo = sub.add_parser(
        "demo",
        help="Executa demonstração rápida (< 10s) de predição e explicabilidade SHAP em um isolado FASTA (Seção 16).",
    )
    p_demo.add_argument("--especie", default="Escherichia coli", help="Espécie bacteriana para demonstração")
    p_demo.add_argument("--output-dir", default="data/processed", help="Diretório dos modelos treinados")

    p_dl = sub.add_parser(
        "download-bvbrc",
        help="Baixa fenótipos de resistência (MIC) da API REST pública do BV-BRC com manifesto SHA-256 (Seção 3.1).",
    )
    p_dl.add_argument("--especie", default="Escherichia coli", help="Espécie bacteriana alvo")
    p_dl.add_argument("--limite", type=int, default=250, help="Limite de registros na consulta")
    p_dl.add_argument("--raw-dir", default="data/raw", help="Diretório para salvar dados brutos")

    for etapa in ("ingest", "qc", "features"):
        p_et = sub.add_parser(etapa, help=f"Executa etapa '{etapa}' do workflow Snakemake.")
        p_et.add_argument("--input", default="", help="Arquivo de entrada opcional")
        p_et.add_argument("--output-dir", default="data/processed", help="Diretório de saída")

    args = parser.parse_args(argv)
    if args.comando == "download-bvbrc":
        regs_dl, cam_dl, hash_dl = baixar_dados_bvbrc_api(
            especie=args.especie,
            limite=args.limite,
            diretorio_raw=Path(args.raw_dir),
        )
        print(
            json.dumps(
                {
                    "especie": args.especie,
                    "total_registros": len(regs_dl),
                    "arquivo_salvo": str(cam_dl),
                    "sha256": hash_dl,
                },
                indent=2,
                ensure_ascii=False,
            )
        )
        return 0

    if args.comando == "demo":
        caminho_pkl = Path(args.output_dir) / "modelos_servico.pkl"
        if caminho_pkl.exists():
            registro = RegistroModelosServico()
            registro.carregar(caminho_pkl)
        else:
            registro, _ = executar_pipeline_completo(
                output_dir=args.output_dir,
                reports_dir="reports",
                n_isolados_base=120,
                n_boot=20,
            )
        pasta_tmp = Path(args.output_dir) / "demo_isolado"
        pasta_tmp.mkdir(parents=True, exist_ok=True)
        fasta_demo = pasta_tmp / "isolado_demo.fasta"
        if args.especie == "Staphylococcus aureus":
            genes_demo = ["gene_mecA", "mut_grlA_S80F", "mut_gyrA_S84L"]
        elif args.especie == "Klebsiella pneumoniae":
            genes_demo = ["gene_blaCTX_M_15", "gene_qnrS1", "gene_aac_3_IIa"]
        else:
            genes_demo = ["mut_gyrA_S83L", "mut_gyrA_D87N", "mut_parC_S80I", "gene_blaCTX_M_15"]
        conteudo_fasta = sintetizar_fasta_isolado(
            genome_id="DEMO_ISOLADO_01",
            atributos_ativos=genes_demo,
        )
        fasta_demo.write_text(conteudo_fasta, encoding="utf-8")
        stats = validar_fasta_entrada(fasta_demo)
        qc_res = avaliar_qc_genoma(
            genome_id="DEMO_ISOLADO_01",
            tamanho_pb=int(stats["tamanho_pb"]),
            contigs=int(stats["contigs"]),
            n50=int(stats["n50"]),
            fracao_n=float(stats["fracao_n"]),
            especie_declarada=args.especie,
            especie_inferida=args.especie,
        )
        predicoes = registro.predizer_isolado_fasta(fasta_demo, especie=args.especie)
        saida_demo = {
            "modo": "demo_interativa",
            "especie": args.especie,
            "arquivo_fasta": str(fasta_demo),
            "qc": qc_res.model_dump(),
            "predicoes": predicoes,
        }
        print(json.dumps(saida_demo, indent=2, ensure_ascii=False))
        return 0

    if args.comando in {"run-all", "ingest", "qc", "features", None}:
        out_dir = getattr(args, "output_dir", "data/processed")
        rep_dir = getattr(args, "reports_dir", "reports")
        n_iso = getattr(args, "n_isolados", 360)
        n_boot = getattr(args, "n_boot", 200)
        _, relatorio = executar_pipeline_completo(
            output_dir=out_dir,
            reports_dir=rep_dir,
            n_isolados_base=n_iso,
            n_boot=n_boot,
        )
        resumo = {
            "status": "sucesso",
            "especie": relatorio["especie_piloto"],
            "total_aprovados_pos_qc_dedup": relatorio["funil_exclusao"]["total_aprovados_pos_qc_dedup"],  # type: ignore[index]
            "relatorio_json": str(Path(rep_dir) / "relatorio_validacao.json"),
            "relatorio_markdown": str(Path(rep_dir) / "relatorio_cientifico.md"),
        }
        print(json.dumps(resumo, indent=2, ensure_ascii=False))
        return 0

    parser.print_help()
    return 1


if __name__ == "__main__":
    sys.exit(main())

