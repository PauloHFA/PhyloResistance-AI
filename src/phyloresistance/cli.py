"""Interface de Linha de Comando (CLI) para reprodução do pipeline PhyloResistance AI."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from phyloresistance.pipeline import executar_pipeline_completo


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="phyloresistance",
        description="PhyloResistance AI — Pipeline reproduzível de predição de AMR e XAI genômico.",
    )
    sub = parser.add_subparsers(dest="comando")

    p_run = sub.add_parser("run-all", help="Executa o pipeline completo ponta a ponta (Fases 1 a 5).")
    p_run.add_argument("--output-dir", default="data/processed", help="Diretório de saída dos artefatos")
    p_run.add_argument("--reports-dir", default="reports", help="Diretório de saída dos relatórios")
    p_run.add_argument("--n-isolados", type=int, default=360, help="Tamanho da coorte piloto")
    p_run.add_argument("--n-boot", type=int, default=200, help="Iterações de bootstrap em blocos")

    for etapa in ("ingest", "qc", "features"):
        p_et = sub.add_parser(etapa, help=f"Executa etapa '{etapa}' do workflow Snakemake.")
        p_et.add_argument("--input", default="", help="Arquivo de entrada opcional")
        p_et.add_argument("--output-dir", default="data/processed", help="Diretório de saída")

    args = parser.parse_args(argv)
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

