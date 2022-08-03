rule star_align_chimera:
    input:
        fq1="fastp/trimmed/{sample}.1.fastq",
        fq2="fastp/trimmed/{sample}.2.fastq",
        idx=config["star"]["index"],
    output:
        chim_junc=temp("star/{sample}/chimera/{sample}.Chimeric.out.junction"),
        bam=temp("star/{sample}/chimera/{sample}.unsorted.bam"),
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
        extra=lambda wildcards: f"--outSAMattrRGline 'ID:{wildcards.sample}\tSM:{wildcards.sample}\tPU:{wildcards.sample}\tPL:ILLUMINA\tCN:IGR\tDS:WES\tPG:STAR' {config['star'].get('chimera_extra')}",
    wrapper:
        "bio/star/align"
