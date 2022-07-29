rule star_align_chimera:
    input:
        fq1="fastp/trimmed/{sample}.1.fastq",
        fq2="fastp/trimmed/{sample}.2.fastq",
        idx=config["star"]["index"],
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
    retries: 2
    log:
        "logs/star/{sample}.log",
    params:
        idx=config["star"]["index"],
        extra=config["star"].get(
            "chimera_extra",
        (
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
                "--outSAMattributes Standard "
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
        ),
    wrapper:
        "bio/star/align"
