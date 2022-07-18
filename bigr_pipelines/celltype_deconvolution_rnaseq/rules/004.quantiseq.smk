rule quantiseq:
    input:
        expr_mat="immunedeconv/TPM.tsv"
    output:
        histogram="quantiseq/celltypes.hist.png",
        dotplot="quantiseq/celltypes.dotplot.png",
        tsv="quantiseq/celltypes.tsv",
        rds="quantiseq/celltypes.RDS",
        plotdir=directory("quantiseq/celltypes.dotplots")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 2048,
        time_min=lambda wildcards, attempt: attempt * 15,
        tmp="tmp"
    message:
        "Using QuantiSeq to deconvolute expression into cell types"
    params:
        gene_col="Hugo_ID"
    log:
        "logs/immunedeconv/quantiseq.log"
    wrapper:
        "bio/immunedeconv/quantiseq"