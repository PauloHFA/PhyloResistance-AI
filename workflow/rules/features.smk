# Regra Snakemake 4: Extração de atributos F1 (genes AMR) e F2 (mutações pontuais) (Seção 5.1)

rule extract_f1_f2_features:
    input:
        clusters="data/processed/clusters_filogeneticos.parquet"
    output:
        f1_f2="data/processed/atributos_f1_f2.parquet"
    conda:
        "../envs/bioinfo.yaml"
    shell:
        "python -m phyloresistance.cli features --output-dir data/processed"

