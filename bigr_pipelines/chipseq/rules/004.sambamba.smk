rule sambamba_markdup:
    input:
        "bowtie2/raw/{sample}.bam"
        "bowtie2/raw/{sample}.bam.bai"
    output:
        temp("sambamba/markdup/{sample}.bam")
    threads: 5
    resources:
        mem_mb=get_2gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir="tmp"
    log:
        "logs/sambamba/{sample}.bwa.log"
    params:
        "--remove-duplicates"
    wrapper:
        "bio/sambamba/markdup"


rule samtools_index_sambamba:
    input:
        "sambamba/markdup/{sample}.bam"
    output:
        temp("sambamba/markdup/{sample}.bam.bai")
    threads: 4
    resources:
        mem_mb=get_2gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir="tmp"
    log:
        "log/samtools/{sample}.index.bwa.log"
    params:
        extra=""
    wrapper:
        "bio/samtools/index"