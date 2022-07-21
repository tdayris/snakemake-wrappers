rule sambamba_markdup:
    input:
        "bowtie2/sorted/{sample}.bam",
        "bowtie2/sorted/{sample}.bam.bai",
    output:
        temp("sambamba/markdup/{sample}.bam"),
    threads: 20
    resources:
        mem_mb=get_2gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir="tmp",
    log:
        "logs/sambamba/{sample}.bwa.log",
    params:
        "--remove-duplicates",
    wrapper:
        "bio/sambamba/markdup"


rule sambamba_index_sambamba:
    input:
        "sambamba/markdup/{sample}.bam",
    output:
        temp("sambamba/markdup/{sample}.bam.bai"),
    threads: 8
    resources:
        mem_mb=get_2gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir="tmp",
    log:
        "logs/sambamba/{sample}.index.bwa.log",
    params:
        extra="",
    wrapper:
        "bio/sambamba/index"
