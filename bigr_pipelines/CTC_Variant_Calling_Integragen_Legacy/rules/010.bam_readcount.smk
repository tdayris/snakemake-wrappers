"""
 bam-readcount v0.8 : bam-readcount -b 20 -q 20 -w 1 -f hg38_basechr.fa {sample}.cutadapt.sorted.bam >{sample}.cutadapt.sorted.bam.readcount
 """
rule bam_readcount:
    input:
        ref = config["ref"]["fasta"],
        bam = "sambamba/sort/{sample}.{status}.bam"
    output:
        temp("bam_readcount/{sample}.{status}.tsv")
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 1024 * 6,
        time_min = lambda wildcards, attempt: attempt * 35,
        tmpdir = "tmp"
    log:
        "logs/bam_readcount/{sample}.{status}.log"
    conda:
        "envs/bam_readcount.yaml"
    params:
        " -b 20 -q 20 -w 1"
    shell:
        " bam-readcount {params} -f {input.ref} {input.bam} > {output} 2> {log}"
