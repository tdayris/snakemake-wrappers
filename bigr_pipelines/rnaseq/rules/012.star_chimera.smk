rule star_align_chimera:
    input:
        fq1="fastp/trimmed/{sample}.1.fastq",
        fq2="fastp/trimmed/{sample}.2.fastq",
        index=config["star"]["index"],
    output:
        chim_junc=temp("star/{sample}/chimera/{sample}.Chimeric.out.junction"),
        bam=temp("star/{sample}/chimera/{sample}.bam"),
        sj=temp("star/{sample}/chimera/{sample}.SJ.out.tab"),
        log=temp("star/{sample}/chimera/{sample}.Log.out"),
        log_progress=temp("star/{sample}/chimera/{sample}.Log.progress.out"),
        log_final=temp("star/{sample}/chimera/{sample}.Log.final.out"),
    threads: config.get("max_threads", 20)
    resources:
        mem_mb=get_75gb_and_5gb_per_attempt,
        time_min=get_90min_per_attempt,
        tmpdir="tmp",
    log:
        "logs/star/{sample}.log",
    params:
        extra=(
            "--outReadsUnmapped None "
            "--twopassMode Basic "
            "--outSAMstrandField intronMotif "
            "--outSAMunmapped Within "
            "--chimSegmentMin 12 "
            "--chimJunctionOverhangMin 8 "
            "--chimOutJunctionFormat 1 "
            "--alignSJDBoverhangMin 10 "
            "--alignMatesGapMax 100000 "
            "--alignIntronMax 100000 "
            "--alignSJstitchMismatchNmax 5 -1 5 5 "
            "--outSAMattrRGline ID:GRPundef "
            "--chimMultimapScoreRange 3 "
            "--chimScoreJunctionNonGTAG -4 "
            "--chimMultimapNmax 20 "
            "--chimNonchimScoreDropMin 10 "
            "--peOverlapNbasesMin 12 "
            "--peOverlapMMp 0.1 "
            "--alignInsertionFlush Right "
            "--alignSplicedMateMapLminOverLmate 0 "
            "--alignSplicedMateMapLmin 30 "
        ),
    wrapper:
        "bio/star/align"


rule samtools_index_chimera_bam:
    input:
        "star/{sample}/chimera/{sample}.bam",
    output:
        temp("star/{sample}/chimera/{sample}.bam.bai"),
    threads: min(config.get("max_threads", 20), 4)
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_15min_per_attempt,
        tmpdir="tmp",
    params:
        extra="",
    log:
        "logs/samtools/index/{sample}.chimera.log",
    wrapper:
        "bio/samtools/index"


rule picard_collect_multiple_metrics_chimera:
    input:
        bam="star/{sample}/chimera/{sample}.bam",
        bai="star/{sample}/chimera/{sample}.bam.bai",
        ref=config["reference"]["genome"],
        ref_idx=config["reference"]["genome_index"],
    output:
        temp(
            multiext(
                "picard/stats/{sample}",
                ".alignment_summary_metrics",
                ".insert_size_metrics",
                ".insert_size_histogram.pdf",
                ".quality_distribution_metrics",
                ".quality_distribution.pdf",
                ".gc_bias.detail_metrics",
                ".gc_bias.summary_metrics",
                ".gc_bias.pdf",
                ".rna_metrics",
                ".bait_bias_detail_metrics",
                ".bait_bias_summary_metrics",
                ".error_summary_metrics",
                ".pre_adapter_detail_metrics",
                ".pre_adapter_summary_metrics",
            )
        ),
    threads: 1
    resources:
        mem_mb=get_4gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir="tmp",
    log:
        "logs/picard/multiple_metrics/{sample}.chimera.log",
    params:
        extra=(
            "--VALIDATION_STRINGENCY LENIENT "
            "--METRIC_ACCUMULATION_LEVEL null "
            "--METRIC_ACCUMULATION_LEVEL SAMPLE "
            "--REF_FLAT ref_flat.txt"
        ),
    wrapper:
        "bio/picard/collectmultiplemetrics"


rule samtools_cram_chimera:
    input:
        "star/{sample}/chimera/{sample}.bam",
        bam_index="star/{sample}/chimera/{sample}.bam.bai",
        ref=config["reference"]["genome"],
        ref_idx=config["reference"]["genome_index"],
    output:
        "star/{sample}/chimera/{sample}.cram",
    threads: min(config.get("max_threads", 2), 2)
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_15min_per_attempt,
        tmpdir="tmp",
    params:
        extra="",
    log:
        "logs/samtools/cram/{sample}.chimera.log",
    wrapper:
        "bio/samtools/view"
