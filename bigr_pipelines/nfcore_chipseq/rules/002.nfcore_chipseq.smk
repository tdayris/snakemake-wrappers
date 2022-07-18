###############################
###   Run nf-core chipseq   ###
###############################


rule ncfore_chipseq_pipeline:
    input:
        input="design.csv",
        fasta=config["ref"]["fasta"],
        gtf=config["ref"]["gtf"],
        c="nextflow.config"
    output:
        "results/mulitqc/broadPeak/multiqc_report.html"
    threads: min(config["threads"], 10)
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 8,
        time_min=lambda wildcards, attempt: attempt * 60 * 8,
        tmpdir="tmp"
    envmodules:
        "java/12.0.2",
        "nextflow/20.11.0-edge",
        "singularity/3.6.3"
    params:
        pipeline="nf-core/chipseq",
        revision="1.2.2",
        profile="singularity",
        genome="GRCh38",
        igenomes_base=config["igenome_base"]
    handover: True
    log:
        "logs/nf-core/chipseq_pipeline.log"
    wrapper:
        "v1.7.0/utils/nextflow"