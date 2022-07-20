"""
This snakefile contains all rules related to immunedeconv
"""


# Call this fantom rule to access all results from immunedeconv
rule immunedeconv_results:
    input:
        graphs=expand(
            "{tool}/celltypes.{graph}.png",
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
    retries: 2
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
        histogram="xcell/celltypes.hist.png",
        dotplot="xcell/celltypes.dotplot.png",
        tsv="xcell/celltypes.tsv",
        rds="xcell/celltypes.RDS",
        plotdir=directory("xcell/celltypes.dotplots"),
    threads: 1
    resources:
        mem_mb=get_4gb_per_attempt,
        time_min=get_15min_per_attempt,
        tmp="tmp",
    retries: 2
    message:
        "Using xCell to deconvolute expression into cell types"
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
        histogram="quantiseq/celltypes.hist.png",
        dotplot="quantiseq/celltypes.dotplot.png",
        tsv="quantiseq/celltypes.tsv",
        rds="quantiseq/celltypes.RDS",
        plotdir=directory("quantiseq/celltypes.dotplots"),
    threads: 1
    resources:
        mem_mb=get_4gb_per_attempt,
        time_min=get_15min_per_attempt,
        tmp="tmp",
    retries: 2
    message:
        "Using QuantiSeq to deconvolute expression into cell types"
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
        histogram="epic/celltypes.hist.png",
        dotplot="epic/celltypes.dotplot.png",
        tsv="epic/celltypes.tsv",
        rds="epic/celltypes.RDS",
        plotdir=directory("epic/celltypes.dotplots"),
    threads: 1
    resources:
        mem_mb=get_4gb_per_attempt,
        time_min=get_15min_per_attempt,
        tmp="tmp",
    retries: 2
    params:
        gene_col="Hugo_ID",
    message:
        "Using EPIC to deconvolute expression into cell types"
    log:
        "logs/immunedeconv/epic.log",
    wrapper:
        "bio/immunedeconv/epic"


# Use mcpcounter to perform immune deconvolution
rule mcpcounter:
    input:
        expr_mat="immunedeconv/TPM.tsv",
    output:
        histogram="mcpcounter/celltypes.hist.png",
        dotplot="mcpcounter/celltypes.dotplot.png",
        tsv="mcpcounter/celltypes.tsv",
        rds="mcpcounter/celltypes.RDS",
        plotdir=directory("mcpcounter/celltypes.dotplots"),
    threads: 1
    resources:
        mem_mb=get_4gb_per_attempt,
        time_min=get_15min_per_attempt,
        tmp="tmp",
    retries: 2
    message:
        "Using MCP-Counter to deconvolute expression into cell types"
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
        histogram="cibersort_abs/celltypes.hist.png",
        dotplot="cibersort_abs/celltypes.dotplot.png",
        tsv="cibersort_abs/celltypes.tsv",
        rds="cibersort_abs/celltypes.RDS",
        plotdir=directory("cibersort_abs/celltypes.dotplots"),
    message:
        "Using Cibersort-absolute to deconvolute expression into cell types"
    threads: 1
    resources:
        mem_mb=get_4gb_per_attempt,
        time_min=get_15min_per_attempt,
        tmp="tmp",
    retries: 2
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
        histogram="cibersort/celltypes.hist.png",
        dotplot="cibersort/celltypes.dotplot.png",
        tsv="cibersort/celltypes.tsv",
        rds="cibersort/celltypes.RDS",
        plotdir=directory("cibersort/celltypes.dotplots"),
    threads: 1
    resources:
        mem_mb=get_4gb_per_attempt,
        time_min=get_15min_per_attempt,
        tmp="tmp",
    retries: 2
    message:
        "Using cibersort to deconvolute expression into cell types"
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
    retries: 2
    message:
        "Gathering Cibersort requirements"
    log:
        "logs/cibersort.bin.log",
    conda:
        "envs/bash.yaml"
    params:
        rsync="--verbose --checksum --human-readable --partial --progress",
        chmod="u+x",
    shell:
        "rsync {params.rsync} {input.cibersort_binary} {output.cibersort_binary} > {log} 2>&1 && "
        "chmod {params.chmod} {output.cibersort_binary} >> {log} 2>&1 && "
        "rsync {params.rsync} {input.cibersort_mat} {output.cibersort_mat} >> {log} 2>&1"
