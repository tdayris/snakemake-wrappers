###############################
### Prepare nf-core chipseq ###
###############################


rule

rule nextflow_config:
    output:
        "nextflow.config"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 512,
        time_min=lambda wildcards, attempt: attempt * 5,
        tmpdir="tmp"
    params:
        text=default_config
    log:
        "logs/nf-core/config.log"
    shell:
        "echo {params.text} > {output} 2> {log}"