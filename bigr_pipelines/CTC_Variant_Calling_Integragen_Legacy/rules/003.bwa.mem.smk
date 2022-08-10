"""
No command provided
"""


rule bwa_index:
    input:
        config["ref"]["fasta"],
    output:
        multiext(
            "../tmp/bwa/index/GRCh38.99.homo_sapiens",
            ".amb",
            ".ann",
            ".bwt",
            ".pac",
            ".sa",
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 1024 * 10 * attempt,
        time_min=lambda wildcards, attempt: 120 * attempt,
        tmpdir="tmp",
    envmodules:
        "bwa/0.7.15",
    log:
        "bwa/index.log",
    params:
        "",
    shell:
        "bwa index {params} {input} -p bwa/index/GRCh38.99.homo_sapiens > {log} 2>&1"


"""
bwa 0.7.15-r1140 : bwa mem -R '@rg\tID:GRCh38\tSM:{sample}\tPL:Illumina' -t {threads} hg38_basechr.fa {sample}.cutadapt.fastq >{sample}.cutadapt.sam
"""


rule bwa_mem:
    input:
        fq="cutadapt/{sample}.{status}.fastq",
        index=config.get(
            "bwa_index",
            multiext(
                "bwa/index/GRCh38.99.homo_sapiens",
                ".amb",
                ".ann",
                ".bwt",
                ".pac",
                ".sa",
            ),
        ),
    output:
        temp("bwa/mem/{sample}.{status}.sam"),
    threads: 20
    resources:
        mem_mb=lambda wildcards, attempt: 1024 * 68 * attempt,
        time_min=lambda wildcards, attempt: 120 * attempt,
        tmpdir="tmp",
    group:
        "map_n_sort"
    envmodules:
        "bwa/0.7.15",
    log:
        "logs/bwa/mem/{sample}.{status}.log",
    params:
        " -R '@rg\tID:GRCh38\tSM:{sample}\tPL:Illumina'",
    shell:
        "bwa mem {params} -t {threads} /bwa/index/GRCh38.99.homo_sapiens {input.fq} > {output} 2>&1"
