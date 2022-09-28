# Estimate read distribution among dna strands
"""
019.rseqc_read_distribution
from:
-> 010.star_align_variants
-> 012.star_align_chimera
by:
-> Snakefile.star_fusion_results
"""


rule rseqc_read_distribution:
    input:
        bam="010.star/{sample}/{maptype}/{sample}.bam",
        bai="010.star/{sample}/{maptype}/{sample}.bam.bai",
        refgene=config["reference"]["refgene_model"],
    output:
        temp("019.rseqc/read_distribution/{maptype}/{sample}.read_distribution.tsv"),
    threads: 1
    resources:
        mem_mb=get_4gb_per_attempt,
        time_min=get_2h_per_attempt,
        tmpdir="tmp",
    retries: 1
    log:
        "logs/019.rseqc/read_distribution/{sample}.{maptype}.log",
    conda:
        str(workflow_source_dir / "envs" / "rseqc.yaml")
    shell:
        "read_distribution.py "
        "--input-file {input.bam} "
        "--refgene {input.refgene} "
        "> {output} 2> {log} "


# Compute transcript integrity number
"""
019.rseqc_tin
from:
-> 010.star_align_variants
-> 012.star_align_chimera
by:
-> Snakefile.star_fusion_results
"""


rule rseqc_tin:
    input:
        bam="010.star/{sample}/{maptype}/{sample}.bam",
        bai="010.star/{sample}/{maptype}/{sample}.bam.bai",
        refgene=config["reference"]["refgene_model"],
    output:
        temp("019.rseqc/tin/{maptype}/{sample}.summary.txt"),
    threads: 1
    resources:
        mem_mb=get_2gb_per_attempt,
        time_min=get_10h_per_attempt,
        tmpdir="tmp",
    retries: 1
    log:
        "logs/019.rseqc/tin/{sample}.{maptype}.log",
    conda:
        str(workflow_source_dir / "envs" / "rseqc.yaml")
    params:
        extra=config["rseqc"].get("tin", "--sample-size 100 --minCov 10"),
    shell:
        "tin.py "
        "--input {input.bam} "
        "--refgene {input.refgene} "
        "{params.extra} "
        "> {output} 2> {log} "


# Acquire multiple stats over BAM file
"""
019.rseqc_bam_stat
from:
-> 010.star_align_variants
-> 012.star_align_chimera
by:
-> Snakefile.star_fusion_results
"""


rule rseqc_bam_stat:
    input:
        bam="010.star/{sample}/{maptype}/{sample}.bam",
        bai="010.star/{sample}/{maptype}/{sample}.bam.bai",
    output:
        temp("019.rseqc/bam_stat/{maptype}/{sample}.txt"),
    threads: 1
    resources:
        mem_mb=get_4gb_per_attempt,
        time_min=get_2h_per_attempt,
        tmpdir="tmp",
    retries: 1
    log:
        "logs/019.rseqc/bam_stat/{sample}.{maptype}.log",
    conda:
        str(workflow_source_dir / "envs" / "rseqc.yaml")
    params:
        extra=config["rseqc"].get("bam_stat", "--mapq 30"),
    shell:
        "bam_stat.py --input-file {input.bam} > {output} 2> {log}"


# Compute RNA coverage bias over gene bodies
"""
019.rseqc_gene_body_coverage
from
-> 010.star_align_variants
-> 012.star_align_chimera
by
-> Snakefile.star_fusion_results
"""


rule seqc_gene_body_coverage:
    input:
        bam="010.star/{sample}/{maptype}/{sample}.bam",
        bai="010.star/{sample}/{maptype}/{sample}.bam.bai",
        refgene=config["reference"]["refgene_model"],
    output:
        txt=temp("019.rseqc/gene_body_coverage/{maptype}/{sample}.geneBodyCoverage.txt"),
        pdf=temp(
            "019.rseqc/gene_body_coverage/{maptype}/{sample}.geneBodyCoverage.curves.pdf"
        ),
        r=temp("019.rseqc/gene_body_coverage/{maptype}/{sample}.geneBodyCoverage.r"),
    threads: 1
    resources:
        mem_mb=get_2gb_per_attempt,
        time_min=get_16h_per_attempt,
        tmpdir="tmp",
    retries: 1
    log:
        "logs/019.rseqc/gene_body_coverage/{sample}.{maptype}.log",
    conda:
        str(workflow_source_dir / "envs" / "rseqc.yaml")
    params:
        extra=config["rseqc"].get("gene_body", "--format pdf"),
        prefix=lambda wildcards: f"rseqc/gene_body_coverage/{wildcards.maptype}/{wildcards.sample}",
    shell:
        "geneBody_coverage.py "
        "--input {input.bam} "
        "--refgene {input.refgene} "
        "--out-prefix {params.prefix} "
        "{params.extra} "
        "> {log} 2>&1 "


# Aggregate stats over junction coverage
"""
019.seqc_junction_annotation
from
-> 010.star_align_variants
-> 012.star_align_chimera
by
-> Snakefile.star_fusion_results
"""


rule rseqc_junction_annotation:
    input:
        bam="010.star/{sample}/{maptype}/{sample}.bam",
        bai="010.star/{sample}/{maptype}/{sample}.bam.bai",
        refgene=config["reference"]["refgene_model"],
    output:
        txt=temp("019.rseqc/junction_annotation/{maptype}/{sample}.txt"),
        pdf_splice=temp(
            "019.seqc/junction_annotation/{maptype}/{sample}.splice_junction.pdf"
        ),
        pdf_events=temp(
            "019.rseqc/junction_annotation/{maptype}/{sample}.splice_events.pdf"
        ),
        xls=temp("019.rseqc/junction_annotation/{maptype}/{sample}.junction.xls"),
        bed=temp("019.rseqc/junction_annotation/{maptype}/{sample}.junction.bed"),
        rscript=temp("rseqc/junction_annotation/{maptype}/{sample}.junction_plot.r"),
        interact=temp(
            "019.rseqc/junction_annotation/{maptype}/{sample}.junction.Interact.bed"
        ),
    threads: 1
    resources:
        mem_mb=get_2gb_per_attempt,
        time_min=get_1h_per_attempt,
        tmpdir="tmp",
    retries: 1
    log:
        "logs/019.rseqc/junction_annotation/{sample}.{maptype}.log",
    conda:
        str(workflow_source_dir / "envs" / "rseqc.yaml")
    params:
        extra=config["rseqc"].get("junction_annotation", "--mapq 30 --min-intron 50"),
        prefix=lambda wildcards: f"rseqc/junction_annotation/{wildcards.maptype}/{wildcards.sample}",
    shell:
        "junction_annotation.py "
        "--input-file {input.bam} "
        "--refgene {input.refgene} "
        "--out-prefix {params.prefix} "
        "{params.extra} "
        "> {output.txt} 2> {log} "
