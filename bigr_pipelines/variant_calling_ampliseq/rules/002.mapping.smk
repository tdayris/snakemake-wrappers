# This rule indexes a given bam file, using samtools index. Most of the time,
# fasta indexes are not explicitely requested by softwares, but they will crash
# if this index is missing.

rule samtools_index_bam:
    input:
        "samtools/filter/{sample}.bam"
    output:
        "samtools/filter/{sample}.bam.bai"
    message:
        "Indexing {wildcards.sample}.bam with Samtools index"
    threads: 20
    resources:
        mem_mb=lambda wildcards, attempt: min(attempt * 1024, 8192),
        time_min=lambda wildcards, attempt: attempt * 30,
        tmpdir="tmp"
    log:
        "samtools/intdex/{sample}/index.log"
    wrapper:
        "bio/samtools/index"


# Filter a bam over the capturekit bed file
rule samtools_filter_bed:
    input:
        aln="sambamba/sort/{sample}.bam",
        fasta=config["ref"]["fasta"],
        fasta_idx=get_fai(config["ref"]["fasta"]),
        fasta_dict=get_dict(config["ref"]["fasta"]),
        bed=config["ref"]["capture_kit_bed"]
    output:
        temp("samtools/filter/{sample}.bam")
    threads: 10
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 2048,
        time_min=lambda wildcards, attempt: attempt * 15,
        tmpdir="tmp"
    params:
        extra="-bh",
        position=""
    log:
        "logs/samtools/filter/{sample}.log"
    wrapper:
        "bio/samtools/view"


# Load and run bwa mapping meta-wrapper

bwa_meta_config = {
    "threads": config["threads"], 
    "genome": config["ref"]["fasta"]
}

module bwa_fixmate_meta:
    snakefile: "../../../meta/bio/bwa_fixmate/test/Snakefile"
    config: bwa_meta_config


use rule * from bwa_fixmate_meta as bwa_fixmate_meta_*


use rule bwa_mem from bwa_fixmate_meta with:
    input:
        reads=expand(
            "fastp/trimmed/pe/{sample}.{stream}.fastq",
            stream=["1", "2"],
            allow_missing=True
        ),
        index=multiext(
            "bwa_mem2/index/genome", ".0123", ".amb", ".ann", ".pac"
        )