rule uncompress_bam:
    input:
        aln="sambamba/markdup/{sample}.bam",
        aln_idx="sambamba/markdup/{sample}.bam.bai",
        ref=config["reference"]["genome"],
        ref_idx=config["reference"]["genome_index"],
    output:
        temp("samtools/view/{sample}.sam"),
    threads: 4
    resources:
        mem_mb=get_2gb_per_attempt,
        time_min=get_2h_per_attempt,
        tmpdir="tmp",
    log:
        "logs/samtools/view/{sample}.filter.log",
    params:
        extra="-q 2 -F 0x04",
    wrapper:
        "bio/samtools/view"


rule extract_fragment_length:
    input:
        "samtools/view/{sample}.sam",
    output:
        temp("awk/fragments/{sample}.duplicated.unsorted.raw.txt"),
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir="tmp",
    log:
        "logs/awk/fragments/{sample}.duplicated.unsorted.raw.log",
    params:
        extra="-F'\\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {print abs($9)}'",
    conda:
        "envs/bash.yaml"
    shell:
        "awk {params.extra} {input} > {output} 2> {log}"


rule sort_fragment_length:
    input:
        "awk/fragments/{sample}.duplicated.unsorted.raw.txt",
    output:
        temp("awk/fragments/{sample}.duplicated.raw.txt"),
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir="tmp",
    log:
        "logs/sort/fragments/{sample}.duplicated.raw.log",
    params:
        extra="",
    conda:
        "envs/bash.yaml"
    shell:
        "sort {params.extra} {input} > {output} 2> {log}"


rule deduplicate_fragment_length:
    input:
        "awk/fragments/{sample}.duplicated.raw.txt",
    output:
        temp("awk/fragments/{sample}.raw.txt"),
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir="tmp",
    log:
        "logs/sort/fragments/{sample}.raw.log",
    params:
        extra="-c",
    conda:
        "envs/bash.yaml"
    shell:
        "uniq {params.extra} {input} > {output} 2> {log}"


rule halve_lengths:
    input:
        "awk/fragments/{sample}.raw.txt",
    output:
        temp("awk/fragments/{sample}.txt"),
