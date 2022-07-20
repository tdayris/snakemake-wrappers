rule star_align_chimera:
    input:
        fq1="fastp/trimmed/{sample}.1.fq.gz",
        fq2="fastp/trimmed/{sample}.2.fq.gz",
        index=config["star"]["index"]
    output:
        chim_junc="star/{sample}/{sample}.Chimeric.out.junction",
        bam="star/{sample}/{sample}.bam",
        sj="star/{sample}/{sample}.SJ.out.tab",
        log="star/{sample}/{sample}.Log.out",
        log_progress="star/{sample}/{sample}.Log.progress.out",
        log_final="star/{sample}/{sample}.Log.final.out"
    threads: 20
    resources:
        mem_mb=get_75gb_and_5gb_per_attempt,
        time_min=get_90min_per_attempt,
        tmpdir="tmp"
    log:
        "logs/star/{sample}.log"
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
        )
    wrapper:
        "bio/star/align"