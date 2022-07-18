rule mcpcounter:
    input:
        expr_mat="immunedeconv/TPM.tsv"
    output:
        histogram="mcpcounter/celltypes.hist.png",
        dotplot="mcpcounter/celltypes.dotplot.png",
        tsv="mcpcounter/celltypes.tsv",
        rds="mcpcounter/celltypes.RDS",
        plotdir=directory("mcpcounter/celltypes.dotplots")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 2048,
        time_min=lambda wildcards, attempt: attempt * 15,
        tmp="tmp"
    message:
        "Using MCP-Counter to deconvolute expression into cell types"
    params:
        gene_col="Hugo_ID"
    log:
        "logs/immunedeconv/mcpcounter.log"
    wrapper:
        "bio/immunedeconv/mcpcounter"