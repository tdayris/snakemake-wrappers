rule bedtools_bamtobed:
    input:
        bam="samtools/view/{sample}.bam",
        bai="samtools/view/{sample}.bam.bai",
    output:
        temp("bedtools/bamtobed/{sample}.bed")
    threads: 1
    resources:
        mem_mb=get_2gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir="tmp"
    log:
        "bedtools/bam2bed/{sample}.log"
    params:
        extra="-bedpe"
    conda:
        "envs/bedtools.yaml"
    shell:
        "bedtools bamtobed -i {input.bam} {params.extra} > {output} 2> {log}"


rule keep_pairs_from_same_chr:
    input:
        "bedtools/bamtobed/{sample}.bed"
    output:
        temp("bedtools/awk/{sample}.bed")
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir="tmp"
    log:
        "awk/bedtools/{sample}.log"
    params:
        extra="'$1==$4 && $6-$2 < 1000 {print $0}'"
    conda:
        "envs/awk.yaml"
    shell:
        "awk {params.extra} {input} > {output} 2> {log}"


rule extract_fragments:
    input:
        "bedtools/awk/{sample}.bed"
    output:
        temp("bedtools/fragments/{sample}.unsorted.bed")
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir="tmp"
    log:
        "cut/bedtools/{sample}.log"
    params:
        extra="-f 1,2,6"
    conda:
        "envs/bash.yaml"
    shell:
        "cut {params.extra} {input} > {output} 2> {log}"
    

rule sort_fragments:
    input:
        "bedtools/fragments/{sample}.unsorted.bed"
    output:
        temp("bedtools/fragments/{sample}.bed")
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir="tmp"
    log:
        "cut/bedtools/{sample}.log"
    params:
        extra="-k1,1 -k2,2n -k3,3n"
    conda:
        "envs/bash.yaml"
    shell:
        "sort {params.extra} {input} > {output} 2> {log}"