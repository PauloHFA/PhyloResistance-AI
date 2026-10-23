# Regra Snakemake 1: Ingestão e curadoria padronizada de rótulos (Seção 3.3 e 4)

rule ingest_and_curate:
    output:
        rotulos="data/processed/curadoria_rotulos.parquet",
        funil="data/processed/funil_exclusao.json"
    conda:
        "../envs/bioinfo.yaml"
    shell:
        "python -m phyloresistance.cli ingest --output-dir data/processed"

