"""
010.star_align_variants
from
-> 002.fastp_clean
by:
-> sambamba_sort_star
"""
rule star_align_variants:
    input:
        fq1="002.fastp/trimmed/{sample}.1.fastq",
        fq2="002.fastp/trimmed/{sample}.2.fastq",
        idx=config["star"]["index"],
    output:
        chim_junc=temp("010.star/{sample}/variants/{sample}.Chimeric.out.junction"),
        bam=temp("010.star/{sample}/variants/{sample}.unsorted.bam"),
        sj=temp("010.star/{sample}/variants/{sample}.SJ.out.tab"),
        log=temp("010.star/{sample}/variants/{sample}.Log.out"),
        log_progress=temp("010.star/{sample}/variants/{sample}.Log.progress.out"),
        log_final=temp("010.star/{sample}/variants/{sample}.Log.final.out"),
    threads: min(config.get("max_threads", 20), 20)
    resources:
        mem_mb=get_75gb_and_5gb_per_attempt,
        time_min=get_90min_per_attempt,
        tmpdir="tmp",
    retries: 1
    log:
        "logs/010.star/{sample}.log",
    params:
        idx=config["star"]["index"],
        extra=lambda wildcards: f"--outSAMattrRGline 'ID:{wildcards.sample}\tSM:{wildcards.sample}\tPU:{wildcards.sample}\tPL:ILLUMINA\tCN:IGR\tDS:WES\tPG:STAR' {config['star'].get('variant_extra')}",
    wrapper:
        "bio/star/align"


"""
010.sambamba_sort_star
from
-> 010.star_align_variants
-> 012.star_align_chimera
by:
-> 010.gatk_split_n_cigar_reads
"""
rule sambamba_sort_star:
    input:
        mapping="010.star/{sample}/{maptype}/{sample}.unsorted.bam",
    output:
        mapping=temp("010.star/{sample}/{maptype}/{sample}.bam"),
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
        "logs/010.star/sort/{sample}.{maptype}.log",
    wrapper:
        "bio/sambamba/sort"


"""
010.gatk_split_n_cigar_reads
from
-> 002.sambamba_sort_star
by:
-> 016.mutect2
"""
rule gatk_split_n_cigar_reads:
    input:
        bam="010.star/{sample}/{maptype}/{sample}.bam",
        bai="010.star/{sample}/{maptype}/{sample}.bam.bai",
        ref=config["reference"]["genome"],
        ref_idx=config["reference"]["genome_index"],
        ref_dict=config["reference"]["genome_dict"],
        bed=config["reference"]["capture_kit_bed"],
    output:
        bam=temp("010.gatk/splitncigarreads/{sample}.bam"),
        bai=temp("010.gatk/splitncigarreads/{sample}.bam.bai"),
    threads: 1
    resources:
        mem_mb=get_4gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir="tmp",
    retries: 1
    log:
        "logs/010.gatk/splitncigarreads/{sample}.log",
    params:
        extra=config["gatk"].get("splitncigarreads_extra", "--create-output-bam-index"),
    conda:
        str(workflow_source_dir / "envs" / "gatk.yaml")
    shell:
        "gatk SplitNCigarReads -R {input.ref} -I {input.bam} -O {output.bam} {extra} > {log} 2>&1"
