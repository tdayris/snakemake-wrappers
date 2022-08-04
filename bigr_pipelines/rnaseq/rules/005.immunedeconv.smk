"""
This snakefile contains all rules related to immunedeconv
"""


# Call this fantom rule to access all results from immunedeconv
rule immunedeconv_results:
    input:
        graphs=expand(
            "data_output/{tool}/celltypes.{graph}.png",
            tool=config["immunedeconv"].get(
                "tool_list",
                [
                    "xcell",
                    "quantiseq",
                    "epic",
                    "mcpcounter",
                    "cibersort",
                    "cibersort_abs",
                ],
            ),
            graph=["histogram", "dotplot"],
        ),
        tpm_table=expand(
            "data_output/Quantification/TPM.{feature}.tsv",
            feature=["transcripts", "genes"],
        ),
        raw_counts="data_output/Quantification/Raw.genes.tsv",
        qc="data_output/MultiQC/ImmuneDeconv.html",


# Subset salmon original TPM counts to fulfill immenedeconv requirements
rule subset_gene_counts:
    input:
        table="salmon/TPM.genes.tsv",
    output:
        table=temp("immunedeconv/TPM.tsv"),
    message:
        "Formatting counts for ImmuneDeconv"
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_15min_per_attempt,
        tmpdir="tmp",
    retries: 1
    log:
        "logs/immunedeconv/filter_gene_counts.log",
    params:
        drop_column=["target_id"],
        set_index="Hugo_ID",
        drop_duplicated_lines=True,
        keep_index=True,
        drop_duplicated_index=True,
        override_previous_index=True,
        filters=[["Hugo_ID", "!=", "Unknown"]],
    wrapper:
        "bio/pandas/filter_table"


# Run xCell to perform immune deconvolution
rule immunedeconv_xcell:
    input:
        expr_mat="immunedeconv/TPM.tsv",
    output:
        histogram="data_output/xcell/celltypes.hist.png",
        dotplot="data_output/xcell/celltypes.dotplot.png",
        tsv="xcell/celltypes.tsv",
        rds="xcell/celltypes.RDS",
        plotdir=directory("data_output/xcell/celltypes.dotplots"),
    threads: 1
    resources:
        mem_mb=get_4gb_per_attempt,
        time_min=get_15min_per_attempt,
        tmp="tmp",
    retries: 1
    params:
        gene_col="Hugo_ID",
    log:
        "logs/immunedeconv/xcell.log",
    wrapper:
        "bio/immunedeconv/xcell"


# Run Quantiseq to perform immune deconvolution
rule immunedeconv_quantiseq:
    input:
        expr_mat="immunedeconv/TPM.tsv",
    output:
        histogram="data_output/quantiseq/celltypes.hist.png",
        dotplot="data_output/quantiseq/celltypes.dotplot.png",
        tsv="quantiseq/celltypes.tsv",
        rds="quantiseq/celltypes.RDS",
        plotdir=directory("data_output/quantiseq/celltypes.dotplots"),
    threads: 1
    resources:
        mem_mb=get_4gb_per_attempt,
        time_min=get_15min_per_attempt,
        tmp="tmp",
    retries: 1
    params:
        gene_col="Hugo_ID",
    log:
        "logs/immunedeconv/quantiseq.log",
    wrapper:
        "bio/immunedeconv/quantiseq"


# Run Epic to perform immune deconvolution
rule immunedeconv_epic:
    input:
        expr_mat="immunedeconv/TPM.tsv",
    output:
        histogram="data_output/epic/celltypes.hist.png",
        dotplot="data_output/epic/celltypes.dotplot.png",
        tsv="epic/celltypes.tsv",
        rds="epic/celltypes.RDS",
        plotdir=directory("data_output/epic/celltypes.dotplots"),
    threads: 1
    resources:
        mem_mb=get_4gb_per_attempt,
        time_min=get_15min_per_attempt,
        tmp="tmp",
    retries: 1
    params:
        gene_col="Hugo_ID",
    log:
        "logs/immunedeconv/epic.log",
    wrapper:
        "bio/immunedeconv/epic"


# Use mcpcounter to perform immune deconvolution
rule mcpcounter:
    input:
        expr_mat="immunedeconv/TPM.tsv",
    output:
        histogram="data_output/mcp_counter/celltypes.hist.png",
        dotplot="data_output/mcp_counter/celltypes.dotplot.png",
        tsv="mcp_counter/celltypes.tsv",
        rds="mcp_counter/celltypes.RDS",
        plotdir=directory("data_output/mcp_counter/celltypes.dotplots"),
    threads: 1
    resources:
        mem_mb=get_4gb_per_attempt,
        time_min=get_15min_per_attempt,
        tmp="tmp",
    retries: 1
    params:
        gene_col="Hugo_ID",
    log:
        "logs/immunedeconv/mcpcounter.log",
    wrapper:
        "bio/immunedeconv/mcpcounter"


rule cibersort_abs:
    input:
        expr_mat="immunedeconv/TPM.tsv",
        cibersort_binary="CIBERSORT.R",
        cibersort_mat="LM22.txt",
    output:
        histogram="data_output/cibersort_abs/celltypes.hist.png",
        dotplot="data_output/cibersort_abs/celltypes.dotplot.png",
        tsv="cibersort_abs/celltypes.tsv",
        rds="cibersort_abs/celltypes.RDS",
        plotdir=directory("data_output/cibersort_abs/celltypes.dotplots"),
    threads: 1
    resources:
        mem_mb=get_4gb_per_attempt,
        time_min=get_15min_per_attempt,
        tmp="tmp",
    retries: 1
    params:
        gene_col="Hugo_ID",
    log:
        "logs/immunedeconv/cibersort_abs.log",
    wrapper:
        "bio/immunedeconv/cibersort-abs"


# Use civersort to perform immune deconvolution
rule cibersort:
    input:
        expr_mat="immunedeconv/TPM.tsv",
        cibersort_binary="CIBERSORT.R",
        cibersort_mat="LM22.txt",
    output:
        histogram="data_output/cibersort/celltypes.hist.png",
        dotplot="data_output/cibersort/celltypes.dotplot.png",
        tsv="cibersort/celltypes.tsv",
        rds="cibersort/celltypes.RDS",
        plotdir=directory("data_output/cibersort/celltypes.dotplots"),
    threads: 1
    resources:
        mem_mb=get_4gb_per_attempt,
        time_min=get_15min_per_attempt,
        tmp="tmp",
    retries: 1
    log:
        "logs/immunedeconv/Cibersort.log",
    params:
        gene_col="Hugo_ID",
    wrapper:
        "bio/immunedeconv/cibersort"


rule get_cibersort:
    input:
        cibersort_binary=config["immunedeconv"]["cibersort_binary"],
        cibersort_mat=config["immunedeconv"]["cibersort_mat"],
    output:
        cibersort_binary=temp("CIBERSORT.R"),
        cibersort_mat=temp("LM22.txt"),
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_15min_per_attempt,
        tmp="tmp",
    retries: 1
    log:
        "logs/cibersort.bin.log",
    conda:
        str(workflow_source_dir / "envs" / "bash.yaml")
    params:
        rsync="--verbose --checksum --human-readable --partial --progress",
        chmod="u+x",
    shell:
        "rsync {params.rsync} {input.cibersort_binary} {output.cibersort_binary} > {log} 2>&1 && "
        "chmod {params.chmod} {output.cibersort_binary} >> {log} 2>&1 && "
        "rsync {params.rsync} {input.cibersort_mat} {output.cibersort_mat} >> {log} 2>&1"


rule multiqc_config_immunedeconv:
    input:
        expand("cibersort/celltypes.tsv", tool=tool_list),
    output:
        yaml=temp("immunedeconv/multiqc_config_mqc.yaml"),
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_15min_per_attempt,
        tmpdir="tmp",
    log:
        "logs/multiqc_config_immunedeconv.log",
    params:
        prefix="",
    conda:
        str(workflow_source_dir / "envs" / "python.yaml")
    script:
        str(workflow_source_dir / "scripts" / "immunedeconv_to_multiqc.py")


rule multiqc:
    input:
        salmon_quant=lambda wildcards: expand(
            "salmon/pseudo_mapping/{sample}/quant.genes.sf",
            sample=samples_per_prefixes[wildcards.comparison],
        ),
        fastp=lambda wildcards: expand(
            "fastp/{ext}/{sample}.fastp.{ext}",
            sample=samples_per_prefixes[wildcards.comparison],
            ext=["html", "json"],
        ),
        fastq_screen=lambda wildcards: expand(
            "fastq_screen/{sample}.{stream}.fastq_screen.{ext}",
            sample=samples_per_prefixes[wildcards.comparison],
            stream=streams,
            ext=["txt", "png"],
        ),
    output:
        "data_output/MultiQC/ImmuneDeconv.html",
    threads: 1
    resources:
        mem_mb=get_2gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir="tmp",
    log:
        "logs/multiqc/immunedeconv.log",
    wrapper:
        "bio/multiqc"
