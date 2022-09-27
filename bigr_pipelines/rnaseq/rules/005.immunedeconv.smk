"""
This snakefile contains all rules related to immunedeconv
"""


# Subset salmon original TPM counts to fulfill immenedeconv requirements
"""
005.subset_gene_counts
from:
-> 004.aggregate_gene_counts
by:
-> 005.immunedeconv_xcell
-> 005.immunedeconv_quantiseq
-> 005.immunedeconv_epic
-> 005.immunedeconv_mcpcounter
-> 005.immunedeconv_cibersort_abs
-> 005.immunedeconv_cibersort
"""
rule 005_subset_gene_counts:
    input:
        table="data_output/Quantification/TPM.genes.tsv",
    output:
        table=temp("005.immunedeconv/TPM.tsv"),
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_15min_per_attempt,
        tmpdir="tmp",
    retries: 1
    log:
        "logs/005.immunedeconv/filter_gene_counts.log",
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
"""
005.immunedeconv_xcell
from:
-> 005.subset_gene_counts
-> 005.multiqc_config_immunedeconv
"""
rule 005_immunedeconv_xcell:
    input:
        expr_mat="005.immunedeconv/TPM.tsv",
    output:
        histogram=protected("data_output/xcell/celltypes.hist.png"),
        dotplot=protected("data_output/xcell/celltypes.dotplot.png"),
        tsv=temp("005.xcell/celltypes.tsv"),
        rds=temp("005.xcell/celltypes.RDS"),
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
        "logs/005.immunedeconv/xcell.log",
    wrapper:
        "bio/immunedeconv/xcell"


# Run Quantiseq to perform immune deconvolution
"""
005.immunedeconv_quantiseq
from:
-> 005.subset_gene_counts
-> 005.multiqc_config_immunedeconv
"""
rule 005_immunedeconv_quantiseq:
    input:
        expr_mat="005.immunedeconv/TPM.tsv",
    output:
        histogram=protected("data_output/quantiseq/celltypes.hist.png"),
        dotplot=protected("data_output/quantiseq/celltypes.dotplot.png"),
        tsv=temp("005.quantiseq/celltypes.tsv"),
        rds=temp("005.quantiseq/celltypes.RDS"),
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
        "logs/005.immunedeconv/quantiseq.log",
    wrapper:
        "bio/immunedeconv/quantiseq"


# Run Epic to perform immune deconvolution
"""
005.immunedeconv_epic
from:
-> 005.subset_gene_counts
-> 005.multiqc_config_immunedeconv
"""
rule 005_immunedeconv_epic:
    input:
        expr_mat="005.immunedeconv/TPM.tsv",
    output:
        histogram=protected("data_output/epic/celltypes.hist.png"),
        dotplot=protected("data_output/epic/celltypes.dotplot.png"),
        tsv=temp("005.epic/celltypes.tsv"),
        rds=temp("005.epic/celltypes.RDS"),
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
        "logs/005.immunedeconv/epic.log",
    wrapper:
        "bio/immunedeconv/epic"


# Use mcpcounter to perform immune deconvolution
"""
005.immunedeconv_mcpcounter
from:
-> 005.subset_gene_counts
-> 005.multiqc_config_immunedeconv
"""
rule 005_immunedeconv_mcpcounter:
    input:
        expr_mat="005.immunedeconv/TPM.tsv",
    output:
        histogram=protected("data_output/mcp_counter/celltypes.hist.png"),
        dotplot=protected("data_output/mcp_counter/celltypes.dotplot.png"),
        tsv=temp("005.mcp_counter/celltypes.tsv"),
        rds=temp("005.mcp_counter/celltypes.RDS"),
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
        "logs/005.immunedeconv/mcpcounter.log",
    wrapper:
        "bio/immunedeconv/mcpcounter"


"""
005.immunedeconv_cibersort_abs
from:
-> 005.subset_gene_counts
-> 005.multiqc_config_immunedeconv
"""
rule 005_immunedeconv_cibersort_abs:
    input:
        expr_mat="005.immunedeconv/TPM.tsv",
        cibersort_binary="CIBERSORT.R",
        cibersort_mat="LM22.txt",
    output:
        histogram=protected("data_output/cibersort_abs/celltypes.hist.png"),
        dotplot=protected("data_output/cibersort_abs/celltypes.dotplot.png"),
        tsv=temp("005.cibersort_abs/celltypes.tsv"),
        rds=temp("005.cibersort_abs/celltypes.RDS"),
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
        "logs/005.immunedeconv/cibersort_abs.log",
    wrapper:
        "bio/immunedeconv/cibersort-abs"


# Use civersort to perform immune deconvolution
"""
005.immunedeconv_cibersort
from:
-> 005.subset_gene_counts
"""
rule 005_immunedeconv_cibersort:
    input:
        expr_mat="005.immunedeconv/TPM.tsv",
        cibersort_binary="CIBERSORT.R",
        cibersort_mat="LM22.txt",
    output:
        histogram=protected("data_output/cibersort/celltypes.hist.png"),
        dotplot=protected("data_output/cibersort/celltypes.dotplot.png"),
        tsv=temp("005.cibersort/celltypes.tsv"),
        rds=temp("005.cibersort/celltypes.RDS"),
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


"""
005.get_cibersort
from:
-> Entry job
by:
-> 005.immunedeconv_cibersort
-> 005.immunedeconv_cibersort_abs
"""
rule 005_get_cibersort:
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
        "logs/005.cibersort.bin.log",
    conda:
        str(workflow_source_dir / "envs" / "bash.yaml")
    params:
        rsync="--verbose --checksum --human-readable --partial --progress",
        chmod="u+x",
    shell:
        "rsync {params.rsync} {input.cibersort_binary} {output.cibersort_binary} > {log} 2>&1 && "
        "chmod {params.chmod} {output.cibersort_binary} >> {log} 2>&1 && "
        "rsync {params.rsync} {input.cibersort_mat} {output.cibersort_mat} >> {log} 2>&1"



rule 005_multiqc_config_immunedeconv:
    input:
        expand("005.cibersort/celltypes.tsv", tool=tool_list),
    output:
        yaml=temp("005.immunedeconv/multiqc_config_mqc.yaml"),
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_15min_per_attempt,
        tmpdir="tmp",
    log:
        "logs/005.multiqc_config_immunedeconv.log",
    params:
        prefix="",
    conda:
        str(workflow_source_dir / "envs" / "python.yaml")
    script:
        str(workflow_source_dir / "scripts" / "005.immunedeconv_to_multiqc.py")
