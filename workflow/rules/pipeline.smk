# Regras do pipeline bioinformático e de ML (Snakemake)

rule ingest_and_curate:
    output:
        rotulos="data/processed/curadoria_rotulos.parquet",
        funil="data/processed/funil_exclusao.json"
    conda:
        "../envs/bioinfo.yaml"
    shell:
        "python -m phyloresistance.cli ingest --output-dir data/processed"

rule qc_and_dedup:
    input:
        rotulos="data/processed/curadoria_rotulos.parquet"
    output:
        qc="data/processed/qc_genomas.parquet",
        dedup="data/processed/relatorio_desduplicacao.json"
    conda:
        "../envs/bioinfo.yaml"
    shell:
        "python -m phyloresistance.cli qc --input {input.rotulos} --output-dir data/processed"

rule phylogeny_and_features:
    input:
        qc="data/processed/qc_genomas.parquet"
    output:
        clusters="data/processed/clusters_filogeneticos.parquet",
        f1_f2="data/processed/atributos_f1_f2.parquet",
        f3="data/processed/atributos_f3_kmers.npz"
    conda:
        "../envs/bioinfo.yaml"
    shell:
        "python -m phyloresistance.cli features --output-dir data/processed"

rule train_and_evaluate:
    input:
        clusters="data/processed/clusters_filogeneticos.parquet",
        f1_f2="data/processed/atributos_f1_f2.parquet",
        f3="data/processed/atributos_f3_kmers.npz"
    output:
        validacao="reports/relatorio_validacao.json",
        relatorio="reports/relatorio_cientifico.md"
    conda:
        "../envs/bioinfo.yaml"
    shell:
        "python -m phyloresistance.cli run-all --output-dir data/processed --reports-dir reports"
