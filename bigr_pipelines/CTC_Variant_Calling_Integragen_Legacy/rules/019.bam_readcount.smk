"""
 bam-readcount v0.8 : bam-readcount -b 20 -q 20 -w 1 -f hg38_basechr.fa {sample}.cutadapt.sorted.bam >{sample}.cutadapt.sorted.bam.readcount
 """


rule bam_readcount_ctc:
    input:
        ref=config["ref"]["fasta"],
        bam="sambamba/sort/{sample}_{version}_{manip}_{nb}.bam",
    output:
        temp("bam_readcount/{sample}_{version}_{manip}_{nb}.tsv"),
    threads: 1
    resources:
        mem_mb=get_8gb_per_attempt,
        time_min=get_35min_per_attempt,
        tmpdir=tmp,
    log:
        "logs/bam_readcount/{sample}_{version}_{manip}_{nb}.log",
    conda:
        str(workflow_source_dir / "envs" / "readcount.yaml")
    params:
        " -b 13 -q 13 -w 1",
    shell:
        " bam-readcount {params} -f {input.ref} {input.bam} > {output} 2> {log}"



rule bam_readcount_wbc:
    input:
        ref=config["ref"]["fasta"],
        bam="sambamba/sort/{sample}_{version}_{manip}.bam",
    output:
        temp("bam_readcount/{sample}_{version}_{manip}.tsv"),
    threads: 1
    resources:
        mem_mb=get_8gb_per_attempt,
        time_min=get_35min_per_attempt,
        tmpdir=tmp,
    log:
        "logs/bam_readcount/{sample}_{version}_{manip}.log",
    conda:
        str(workflow_source_dir / "envs" / "readcount.yaml")
    params:
        " -b 13 -q 13 -w 1",
    shell:
        " bam-readcount {params} -f {input.ref} {input.bam} > {output} 2> {log}"



rule bam_readcount:
    input:
        ref=config["ref"]["fasta"],
        bam="sambamba/sort/{sample}.baseline.bam",
    output:
        temp("bam_readcount/{sample}.baseline.tsv"),
    threads: 1
    resources:
        mem_mb=get_8gb_per_attempt,
        time_min=get_35min_per_attempt,
        tmpdir=tmp,
    log:
        "logs/bam_readcount/{sample}.baseline.log",
    conda:
        str(workflow_source_dir / "envs" / "readcount.yaml")
    params:
        " -b 13 -q 13 -w 1",
    shell:
        " bam-readcount {params} -f {input.ref} {input.bam} > {output} 2> {log}"