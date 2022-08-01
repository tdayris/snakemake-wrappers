rule star_align_variants:
    input:
        fq1="fastp/trimmed/{sample}.1.fastq",
        fq2="fastp/trimmed/{sample}.2.fastq",
        idx=config["star"]["index"],
    output:
        chim_junc=temp("star/{sample}/variants/{sample}.Chimeric.out.junction"),
        sam=temp("star/{sample}/variants/{sample}.sam"),
        sj=temp("star/{sample}/variants/{sample}.SJ.out.tab"),
        log=temp("star/{sample}/variants/{sample}.Log.out"),
        log_progress=temp("star/{sample}/variants/{sample}.Log.progress.out"),
        log_final=temp("star/{sample}/variants/{sample}.Log.final.out"),
    threads: min(config.get("max_threads", 20), 20)
    resources:
        mem_mb=get_75gb_and_5gb_per_attempt,
        time_min=get_90min_per_attempt,
        tmpdir="tmp",
    retries: 2
    log:
        "logs/star/{sample}.log",
    params:
        idx=config["star"]["index"],
        extra=config["star"].get(
            "variant_extra",
        (
        "--outFilterType BySJout "
                "--outFilterMultimapNmax 20 "
                "--alignSJoverhangMin 8 "
                "--alignSJDBoverhangMin 1 "
                "--outFilterMismatchNmax 999 "
                "--outFilterMismatchNoverReadLmax 0.04 "
                "--alignIntronMin 20 "
                "--alignIntronMax 1000000 "
                "--alignMatesGapMax 1000000 "
                "--outSAMattributes Standard "
                "--twopassMode Basic "
            ),
        ),
    wrapper:
        "bio/star/align"



rule samtools_view_star:
    input:
        "star/{sample}/{maptype}/{sample}.sam",
        fasta=config["reference"]["genome"],
        fasta_idx=config["reference"]["genome_index"],
        fasta_dict=config["reference"]["genome_dict"],
        bed=config["reference"]["capture_kit_bed"],
    output:
        temp("samtools/filter/{sample}_{maptype}.unsorted.bam"),
    threads: min(config.get("max_threads", 20), 10)
    resources:
        mem_mb=get_2gb_per_attempt,
        time_min=get_35min_per_attempt,
        tmpdir="tmp",
    params:
        extra="-h",
    log:
        "logs/sambamba/view/{sample}.{maptype}.raw_star.log"
    wrapper:
        "bio/samtools/view"



rule sambamba_sort_star:
    input:
        mapping="star/{sample}/{maptype}/{sample}.unsorted.bam"
    output:
        mapping=temp("star/{sample}/{maptype}/{sample}.bam")
    threads: min(config.get("max_threads", 20), 10)
    resources:
        mem_mb=get_2gb_per_attempt,
        time_min=get_90min_per_attempt,
        tmpdir="tmp",
    shadow:
        "shallow"
    params:
        extra=config["sambamba"].get("sort_extra", ""),
    log:
        "logs/star/sort/{sample}.{maptype}.log",
    wrapper:
        "bio/sambamba/sort"