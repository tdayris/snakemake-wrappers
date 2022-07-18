rule epic:
    input:
        expr_mat="immunedeconv/TPM.tsv"
    output:
        histogram="epic/celltypes.hist.png",
        dotplot="epic/celltypes.dotplot.png",
        tsv="epic/celltypes.tsv",
        rds="epic/celltypes.RDS",
        plotdir=directory("epic/celltypes.dotplots")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 2048,
        time_min=lambda wildcards, attempt: attempt * 15,
        tmp="tmp"
    params:
        gene_col="Hugo_ID"
    message:
        "Using EPIC to deconvolute expression into cell types"
    log:
        "logs/immunedeconv/epic.log"
    wrapper:
        "bio/immunedeconv/epic"