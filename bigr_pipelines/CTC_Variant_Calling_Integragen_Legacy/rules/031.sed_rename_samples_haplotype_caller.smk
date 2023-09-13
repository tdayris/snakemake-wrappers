"""
No command line provided
"""

rule sed_rename_sample_wbc_hc:
    input:
        "vep/hc/{sample}_{version}_{manip}.orig.tsv"
    output:
        temp("vep/hc/{sample}_{version}_{manip}.renamed.tsv")
    threads: 1
    resources:
        mem_mb=512,
        time_min=10,
        tmpdir=tmp,
    log:
        "logs/vep/rename/hc/{sample}_{version}_{manip}.wbc.log"
    conda:
        "../envs/bash.yaml"
    params:
        replace=(
            "s|{sample}_{version}_{manip}"
            "|{sample}_{version}_{manip}_WBC|g"
        )
    shell:
        "sed '{params.replace}' {input} > {output} 2> {log}"


rule sed_rename_sample_baseline_hc:
    input:
        "vep/hc/{sample}.baseline.orig.tsv"
    output:
        temp("vep/hc/{sample}.baseline.renamed.tsv")
    threads: 1
    resources:
        mem_mb=512,
        time_min=10,
        tmpdir=tmp,
    log:
        "logs/vep/rename/hc/{sample}.baseline.log"
    params:
        replace=(
            "s|{sample}_baseline|{sample}_Germline_DNA|g;s|{sample}.baseline|{sample}_Germline_DNA|g"
        )
    conda:
        "../envs/bash.yaml"
    shell:
        "sed '{params.replace}' {input} > {output} 2> {log}"