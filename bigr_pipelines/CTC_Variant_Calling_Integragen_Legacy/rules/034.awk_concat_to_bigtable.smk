"""
No command line provided
"""

rule concat_to_bigtable:
    input:
        # Mutect2
        expand(
            "vep/mutect/{ctc_sample}.orig.tsv",
            ctc_sample=ctc_list
        ),
        expand(
            "vep/mutect/{wbc_sample}.renamed.tsv",
            wbc_sample=wbc_list
        ),
        # Haplotype Caller
        expand(
            "vep/hc/{ctc_sample}.orig.tsv",
            ctc_sample=ctc_list
        ),
        expand(
            "vep/hc/{wbc_sample}.renamed.tsv",
            wbc_sample=wbc_list
        ),
        expand(
            "vep/hc/{baseline_sample}.baseline.renamed.tsv",
            baseline_sample=baseline_list,
        ),
        # Bam ReadCount
        expand(
            "vep/brc/{ctc_sample}.orig.tsv",
            ctc_sample=ctc_list
        ),
        expand(
            "vep/brc/{wbc_sample}.renamed.tsv",
            wbc_sample=wbc_list
        ),
        expand(
            "vep/brc/{baseline_sample}.baseline.renamed.tsv",
            baseline_sample=baseline_list,
        ),
    output:
        temp("bigtable/raw.tsv"),
    threads: 1
    resources:
        mem_mb=get_10gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir=tmp,
    log:
        "logs/bigtable/raw.log",
    conda:
        str(workflow_source_dir / "envs" / "bash.yaml")
    params:
        "'(NR == 1) || (FNR > 1)'",
    shell:
        "awk {params} {input} > {output} 2> {log}"