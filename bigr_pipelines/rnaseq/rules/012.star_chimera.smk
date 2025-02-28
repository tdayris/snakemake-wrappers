# Align reads with parameters dedicated to chimera and fusions findings
"""
012.star_align_chimera
from
-> 002.fastp_clean
by:
-> 010.sambamba_sort_star
-> 019.rseqc_read_distribution
-> 019.rseqc_tin
-> 019.rseqc_bam_stat
-> 019.rseqc_gene_body_coverage
-> 019.seqc_junction_annotation
"""


rule star_align_chimera:
    input:
        fq1="002.fastp/trimmed/{sample}.1.fastq",
        fq2="002.fastp/trimmed/{sample}.2.fastq",
        idx=config["star"]["index"],
    output:
        chim_junc=temp("010.star/{sample}/chimera/{sample}.Chimeric.out.junction"),
        bam=temp("010.star/{sample}/chimera/{sample}.unsorted.bam"),
        sj=temp("010.star/{sample}/chimera/{sample}.SJ.out.tab"),
        log=temp("010.star/{sample}/chimera/{sample}.Log.out"),
        log_progress=temp("010.star/{sample}/chimera/{sample}.Log.progress.out"),
        log_final=temp("010.star/{sample}/chimera/{sample}.Log.final.out"),
    threads: config.get("max_threads", 20)
    resources:
        mem_mb=get_75gb_and_5gb_per_attempt,
        time_min=get_90min_per_attempt,
        tmpdir="tmp",
    retries: 1
    log:
        "logs/010.star/{sample}.log",
    params:
        idx=config["star"]["index"],
        extra=lambda wildcards: f"--outSAMattrRGline 'ID:{wildcards.sample}\tSM:{wildcards.sample}\tPU:{wildcards.sample}\tPL:ILLUMINA\tCN:IGR\tDS:WES\tPG:STAR' {config['star'].get('chimera_extra')}",
    wrapper:
        "bio/star/align"
