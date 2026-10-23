# Regra Snakemake 5: Contagem de k-mers canônicos (k=31) KMC3 / unitigs esparsos F3 (Seção 5.2)

rule extract_f3_kmers:
    input:
        f1_f2="data/processed/atributos_f1_f2.parquet"
    output:
        f3="data/processed/atributos_f3_kmers.npz"
    conda:
        "../envs/bioinfo.yaml"
    shell:
        "python -m phyloresistance.cli features --output-dir data/processed"

