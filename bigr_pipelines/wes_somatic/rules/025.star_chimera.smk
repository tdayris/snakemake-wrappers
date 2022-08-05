rule star_align_chimera:
    input:
        fq1="fastp/trimmed/{sample}_{status}.1.fastq",
        fq2="fastp/trimmed/{sample}_{status}.2.fastq",
        idx=config["star"]["index"],
    output:
        chim_junc=temp("star/{sample}/chimera/{sample}_{status}.Chimeric.out.junction"),
        bam=temp("star/{sample}/{sample}_{status}.unsorted.bam"),
        sj=temp("star/{sample}/chimera/{sample}_{status}.SJ.out.tab"),
        log=temp("star/{sample}/chimera/{sample}_{status}.Log.out"),
        log_progress=temp("star/{sample}/chimera/{sample}_{status}.Log.progress.out"),
        log_final=temp("star/{sample}/chimera/{sample}_{status}.Log.final.out"),
    threads: config.get("max_threads", 20)
    resources:
        mem_mb=get_75gb_and_5gb_per_attempt,
        time_min=get_90min_per_attempt,
        tmpdir="tmp",
    retries: 1
    log:
        "logs/star/{sample}_{status}.log",
    params:
        idx=config["star"]["index"],
        extra=lambda wildcards: f"--outSAMattrRGline 'ID:{wildcards.sample}\tSM:{wildcards.sample}\tPU:{wildcards.sample}\tPL:ILLUMINA\tCN:IGR\tDS:WES\tPG:STAR' {config['star'].get('chimera_extra')}",
    wrapper:
        "bio/star/align"



rule sambamba_sort_star:
    input:
        mapping="star/{sample}/{sample}_{status}.unsorted.bam",
    output:
        mapping=temp("star/chimera/{sample}_{status}.bam"),
    threads: min(config.get("max_threads", 20), 10)
    resources:
        mem_mb=get_6gb_per_attempt,
        time_min=get_90min_per_attempt,
        tmpdir="tmp",
    retries: 1
    shadow:
        "shallow"
    params:
        extra=config["sambamba"].get("sort_extra", ""),
    log:
        "logs/star/sort/{sample}_{status}.log",
    wrapper:
        "bio/sambamba/sort"