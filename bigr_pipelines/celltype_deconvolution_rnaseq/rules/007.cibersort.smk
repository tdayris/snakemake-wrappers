
rule cibersort_abs:
    input:
        expr_mat="immunedeconv/TPM.tsv",
        cibersort_binary="CIBERSORT.R",
        cibersort_mat="LM22.txt"
    output:
        histogram="cibersort_abs/celltypes.hist.png",
        dotplot="cibersort_abs/celltypes.dotplot.png",
        tsv="cibersort_abs/celltypes.tsv",
        rds="cibersort_abs/celltypes.RDS",
        plotdir=directory("cibersort_abs/celltypes.dotplots")
    message:
        "Using Cibersort-absolute to deconvolute expression into cell types"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 2048,
        time_min=lambda wildcards, attempt: attempt * 15,
        tmp="tmp"
    params:
        gene_col="Hugo_ID"
    log:
        "logs/immunedeconv/cibersort_abs.log"
    wrapper:
        "bio/immunedeconv/cibersort-abs"


rule cibersort:
    input:
        expr_mat="immunedeconv/TPM.tsv",
        cibersort_binary="CIBERSORT.R",
        cibersort_mat="LM22.txt"
    output:
        histogram="cibersort/celltypes.hist.png",
        dotplot="cibersort/celltypes.dotplot.png",
        tsv="cibersort/celltypes.tsv",
        rds="cibersort/celltypes.RDS",
        plotdir=directory("cibersort/celltypes.dotplots")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 2048,
        time_min=lambda wildcards, attempt: attempt * 15,
        tmp="tmp"
    message:
        "Using cibersort to deconvolute expression into cell types"
    log:
        "logs/immunedeconv/Cibersort.log"
    params:
        gene_col="Hugo_ID"
    wrapper:
        "bio/immunedeconv/cibersort"


rule get_cibersort:
    input:
        cibersort_binary="/mnt/beegfs/software/cibersort/1.0.6/CIBERSORT.R",
        cibersort_mat="/mnt/beegfs/software/cibersort/1.0.6/LM22.txt"
    output:
        cibersort_binary=temp("CIBERSORT.R"),
        cibersort_mat=temp("LM22.txt")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
        time_min=lambda wildcards, attempt: attempt * 15,
        tmpdir="tmp"
    message:
        "Gathering Cibersort requirements"
    log:
        "logs/cibersort.bin.log"
    params:
        rsync = "--verbose --checksum --human-readable",
        chmod = "u+x"
    shell:
        "rsync {params.rsync} {input.cibersort_binary} {output.cibersort_binary} > {log} 2>&1 && "
        "chmod {params.chmod} {output.cibersort_binary} >> {log} 2>&1 && "
        "rsync {params.rsync} {input.cibersort_mat} {output.cibersort_mat} >> {log} 2>&1"