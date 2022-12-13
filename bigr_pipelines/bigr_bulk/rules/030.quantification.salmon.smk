rule salmon_quantification:
    input:
        **unpack(get_salmon_index),
        **unpack(get_salmon_fastq),
    output:
        quant="030.quantification/salmon_{genome_build}.{genome_release}/{sample}/quant.sf"
    threads: 20
    resources:
        mem_mb=get_10gb_and_10gb_per_attempt,
        time_min=get_45min_and_20min_per_attempt,
        tmpdir="tmp"
    log:
        "logs/030.quantification.salmon/{sample}.{genome_build}.{genome_release}.smk"
    params:
        extra="--numBootstraps 100 --gcBias --seqBias --posBias"
    wrapper:
        "v1.20.0/bio/salmon/quant"



            salmon_logs=snakemake.input["salmon_logs"],
            salmon_metas=snakemake.input["salmon_meta"],
            salmon_mapping_plot=snakemake.output["mapping_status_mqc"],
            salmon_assigned_fragments=snakemake.output["mapping_counts_mqc"],
            salmon_mapping_rates_general_stat=snakemake.output["salmon_general_table"],

rule salmon_to_multiqc:
    input:
        salmon_logs=expand(""),
        salmon_meta=expand(""),
    output:
        mapping_status_mqc=temp(""),
        mapping_counts_mqc=temp(""),
        salmon_general_table=temp(""),
    