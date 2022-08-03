rule star_align_variants:
    input:
        fq1="fastp/trimmed/{sample}.1.fastq",
        fq2="fastp/trimmed/{sample}.2.fastq",
        idx=config["star"]["index"],
    output:
        chim_junc=temp("star/{sample}/variants/{sample}.Chimeric.out.junction"),
        bam=temp("star/{sample}/variants/{sample}.unsorted.bam"),
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
        extra=lambda wildcards: f"--outSAMattrRGline 'ID:{wildcards.sample}\tSM:{wildcards.sample}\tPU:{wildcards.sample}\tPL:ILLUMINA\tCN:IGR\tDS:WES\tPG:STAR' {config['star'].get('variant_extra')}",
    wrapper:
        "bio/star/align"


rule sambamba_sort_star:
    input:
        mapping="star/{sample}/{maptype}/{sample}.unsorted.bam",
    output:
        mapping=temp("star/{sample}/{maptype}/{sample}.bam"),
    threads: min(config.get("max_threads", 20), 10)
    resources:
        mem_mb=get_2gb_per_attempt,
        time_min=get_90min_per_attempt,
        tmpdir="tmp",
    retries: 2
    shadow:
        "shallow"
    params:
        extra=config["sambamba"].get("sort_extra", ""),
    log:
        "logs/star/sort/{sample}.{maptype}.log",
    wrapper:
        "bio/sambamba/sort"


rule gatk_split_n_cigar_reads:
    input:
        bam="star/{sample}/{maptype}/{sample}.bam",
        bai="star/{sample}/{maptype}/{sample}.bam.bai",
        ref=config["reference"]["genome"],
        ref_idx=config["reference"]["genome_index"],
        ref_dict=config["reference"]["genome_dict"],
        bed=config["reference"]["capture_kit_bed"],
    output:
        bam=temp("gatk/splitncigarreads/{sample}.bam"),
        bai=temp("gatk/splitncigarreads/{sample}.bam.bai"),
    threads: 1
    resources:
        mem_mb=get_4gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir="tmp",
    retries: 1
    log:
        "logs/gatk/splitncigarreads/{sample}.log",
    params:
        extra=config["gatk"].get("splitncigarreads_extra", "--create-output-bam-index"),
    conda:
        str(workflow_source_dir / "envs" / "gatk.yaml")
    shell:
        "gatk SplitNCigarReads -R {input.ref} -I {input.bam} -O {output.bam} {extra} > {log} 2>&1"
