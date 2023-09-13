"""
No command line provided
"""

rule sed_rename_sample_wbc_mutect2:
    input:
        "vep/mutect/{sample}_{version}_{manip}.orig.tsv"
    output:
        temp("vep/mutect/{sample}_{version}_{manip}.renamed.tsv")
    threads: 1
    resources:
        mem_mb=512,
        time_min=10,
        tmpdir=tmp,
    log:
        "logs/vep/rename/mutect/{sample}_{version}_{manip}.wbc.log"
    conda:
        "../envs/bash.yaml"
    params:
        replace=(
            "s|{sample}_{version}_{manip}"
            "|{sample}_{version}_{manip}_WBC|g"
        )
    shell:
        "sed '{params.replace}' {input} > {output} 2> {log}"