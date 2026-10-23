# Regra Snakemake 2: Controle de Qualidade (QUAST / CheckM2 / confirmação de espécie — Seção 3.4)

rule qc_genomes:
    input:
        rotulos="data/processed/curadoria_rotulos.parquet"
    output:
        qc="data/processed/qc_genomas.parquet"
    conda:
        "../envs/bioinfo.yaml"
    shell:
        "python -m phyloresistance.cli qc --input {input.rotulos} --output-dir data/processed"

