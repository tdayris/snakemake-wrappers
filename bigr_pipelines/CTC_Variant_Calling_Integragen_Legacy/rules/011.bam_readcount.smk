"""
 bam-readcount v0.8 : bam-readcount -b 20 -q 20 -w 1 -f hg38_basechr.fa {sample}.cutadapt.sorted.bam >{sample}.cutadapt.sorted.bam.readcount
 """


rule bam_readcount:
    input:
        ref=config["ref"]["fasta"],
        bam="sambamba/sort/{sample}.{status}.bam",
    output:
        temp("bam_readcount/{sample}.{status}.tsv"),
    threads: 1
    resources:
        mem_mb=get_8gb_per_attempt,
        time_min=get_35min_per_attempt,
        tmpdir="tmp",
    log:
        "logs/bam_readcount/{sample}.{status}.log",
    conda:
        str(workflow_source_dir / "envs" / "readcount.yaml")
    params:
        " -b 20 -q 20 -w 1",
    shell:
        " bam-readcount {params} -f {input.ref} {input.bam} > {output} 2> {log}"
