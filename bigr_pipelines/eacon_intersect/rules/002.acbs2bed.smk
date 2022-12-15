rule cut_acbs:
    input:
        acbs="data_input/{sample}.Cut.cbs"
    output:
        bedlike=pipe("bed/{sample}.unsorted.bedlike")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 768,
        time_min=lambda wildcards, attempt: attempt * 10,
        tmpdir="tmp"
    log:
        "logs/acbs/{sample}.cut.log"
    params:
        extra="-f2,3,4,6 -d $'\t'"
    shell:
        "cut {params.extra} {input.acbs} > {output.bedlike} 2> {log}"


rule sed_acbs:
    input:
        bedlike="bed/{sample}.unsorted.bedlike"
    output:
        bed=pipe("bed/{sample}.unsorted.bed")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 768,
        time_min=lambda wildcards, attempt: attempt * 10,
        tmpdir="tmp"
    log:
        "logs/acbs/{sample}.sed.log"
    params:
        extra="'1d'"
    shell:
        "sed {params.extra} {input.acbs} > {output.bed} 2> {log}"


rule sort_acbs:
    input:
        bedlike="bed/{sample}.unsorted.bed"
    output:
        bed=temp("bed/{sample}.bed")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 768,
        time_min=lambda wildcards, attempt: attempt * 10,
        tmpdir="tmp"
    log:
        "logs/acbs/{sample}.sort.log"
    params:
        extra="-k1,1 -k2,2 -k3,3"
    shell:
        "sed {params.extra} {input.acbs} > {output.bed} 2> {log}"
