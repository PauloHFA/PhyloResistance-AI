# Regra Snakemake 3: Desduplicação genômica Mash/skani e Clustering Filogenético (REG-01, REG-02)

rule dedup_and_phylogeny:
    input:
        qc="data/processed/qc_genomas.parquet"
    output:
        dedup="data/processed/relatorio_desduplicacao.json",
        clusters="data/processed/clusters_filogeneticos.parquet"
    conda:
        "../envs/bioinfo.yaml"
    shell:
        "python -m phyloresistance.cli qc --input {input.qc} --output-dir data/processed"

