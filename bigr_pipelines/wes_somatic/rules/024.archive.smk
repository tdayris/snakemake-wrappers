rule archive_list:
    input:
        expand("data_output/archive/VCF/{sample}.bcf", sample=sample_list),
        expand("data_output/archive/TSV/{sample}.tsv.gz", sample=sample_list),
        expand("data_output/archive/CNV/{sample}.tsv.gz", sample=sample_list),
        expand("data_output/archive/{content}.tsv", content=["TMB", "MSI"]),
        "data_output/archive/Somatic_Variant_Calling.html.gz",
        expand(
            "data_output/archive/{sample}_{status}.cram",
            sample=sample_list,
            status=status_list,
        ),


rule bcftools_archive:
    input:
        "data_output/VCF/{sample}.vcf.gz",
        "data_output/VCF/{sample}.vcf.gz.tbi",
    output:
        "data_output/archive/VCF/{sample}.bcf",
    threads: 3
    resources:
        mem_mb=get_2gb_per_attempt,
        time_min=get_5h_per_attempt,
        tmpdir=tmp,
    log:
        "logs/archive/vcf/{sample}.log",
    params:
        extra="--compression-level 9",
    wrapper:
        str(wrapper_prefix / "bio" / "bcftools" / "view")


rule gzip_tsv:
    input:
        "data_output/TSV/{sample}.tsv",
    output:
        "data_output/archive/TSV/{sample}.tsv.gz",
    threads: 1
    resources:
        mem_mb=get_2gb_per_attempt,
        time_min=get_5h_per_attempt,
        tmpdir=tmp,
    conda:
        str(workflow_source_dir / "envs" / "bash.yaml")
    log:
        "logs/archive/tsv/{sample}.log",
    params:
        "--stdout --force --verbose --best",
    shell:
        "gzip {params} {input} > {output} 2> {log"


use rule gzip_tsv as gzip_cnv with:
    input:
        "data_output/CNV/{sample}.tsv",
    output:
        "data_output/archive/CNV/{sample}.tsv.gz",
    log:
        "logs/archive/cnv/{sample}.log",


use rule gzip_tsv as gzip_tmb with:
    input:
        "data_output/TMB.tsv",
    output:
        "data_output/archive/TMB.tsv",
    log:
        "logs/archive/tmb.log",


use rule gzip_tsv as gzip_msi with:
    input:
        "data_output/MSI.tsv",
    output:
        "data_output/archive/MSI.tsv",
    log:
        "logs/archive/msi.log",


use rule gzip_tsv as gzip_html with:
    input:
        "data_output/MultiQC/Somatic_Variant_Calling.html",
    output:
        "data_output/archive/Somatic_Variant_Calling.html.gz",
    log:
        "logs/archive/multiqc.log",


rule cram_mapping:
    input:
        aln="data_output/BAM/{sample}_{status}.bam",
        aln_idx="data_output/BAM/{sample}_{status}.bam.bai",
        fasta=config["reference"]["fasta"],
        fasta_idx=config["reference"]["fasta_index"],
        bed=config["reference"]["capture_kit_bed"],
    output:
        "data_output/archive/{sample}_{status}.cram",
    threads: min(config.get("max_threads", 4), 4)
    resources:
        mem_mb=get_2gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir=tmp,
    log:
        "logs/archive/cram/{sample}.{status}.log",
    params:
        extra="-h",
    wrapper:
        str(wrapper_prefix / "bio" / "samtools" / "view")
