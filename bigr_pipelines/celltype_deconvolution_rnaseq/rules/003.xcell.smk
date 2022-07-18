rule xcell:
    input:
        expr_mat="immunedeconv/TPM.tsv"
    output:
        histogram="xcell/celltypes.hist.png",
        dotplot="xcell/celltypes.dotplot.png",
        tsv="xcell/celltypes.tsv",
        rds="xcell/celltypes.RDS",
        plotdir=directory("xcell/celltypes.dotplots")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 2048,
        time_min=lambda wildcards, attempt: attempt * 15,
        tmp="tmp"
    message:
        "Using xCell to deconvolute expression into cell types"
    params:
        gene_col="Hugo_ID"
    log:
        "logs/immunedeconv/xcell.log"
    wrapper:
        "bio/immunedeconv/xcell"