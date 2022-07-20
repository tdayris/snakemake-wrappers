rule samtools_view_filter:
    input:
        aln="sambamba/markdup/{sample}.bam",
        aln_idx="sambamba/markdup/{sample}.bam.bai",
        ref=config["reference"]["genome"],
        ref_idx=config["reference"]["genome_index"],
    output:
        protected("samtools/view/{sample}.bam")
    threads: 4
    resources:
        mem_mb=get_2gb_per_attempt,
        time_min=get_2h_per_attempt,
        tmpdir="tmp"
    log:
        "logs/samtools/view/{sample}.filter.log"
    params:
        extra="-q 2 -F 0x04"
    wrapper:
        "bio/samtools/view"


rule samtools_index_filtered:
    input:
        "samtools/view/{sample}.bam"
    output:
        protected("samtools/view/{sample}.bam.bai")
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