rule star_align_variants:
    input:
        fq1="fastp/trimmed/{sample}.1.fastq",
        fq2="fastp/trimmed/{sample}.2.fastq",
        idx=config["star"]["index"],
    output:
        chim_junc=temp("star/{sample}/variants/{sample}.Chimeric.out.junction"),
        bam=temp("star/{sample}/variants/{sample}.bam"),
        sj=temp("star/{sample}/variants/{sample}.SJ.out.tab"),
        log=temp("star/{sample}/variants/{sample}.Log.out"),
        log_progress=temp("star/{sample}/variants/{sample}.Log.progress.out"),
        log_final=temp("star/{sample}/variants/{sample}.Log.final.out"),
        reads_per_gene=temp("star/{sample}/variants/{sample}.counts.tsv")
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
                "--quantMode GeneCounts "
                "--outFilterMultimapNmax 20 "
                "--alignSJoverhangMin 8 "
                "--alignSJDBoverhangMin 1 "
                "--outFilterMismatchNmax 999 "
                "--outFilterMismatchNoverReadLmax 0.04 "
                "--alignIntronMin 20 "
                "--alignIntronMax 1000000 "
                "--alignMatesGapMax 1000000 "
                "--outSAMattributes All "
                "--twopassMode Basic "
            ),
        ),
    wrapper:
        "bio/star/align"
