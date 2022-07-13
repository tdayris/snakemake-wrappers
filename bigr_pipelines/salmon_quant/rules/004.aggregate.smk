# This snakefile contains rules to aggregate
# individual counts (aka. per-sample estimates)
# to a single count table.


# Work on raw counts, usefull for manual
# DESeq2/EdgeR/...
rule aggregate_raw_counts:
    input:
        quant=expand(
            "salmon/pseudo_mapping/{sample}/quant.genes.sf", sample=design["Sample_id"]
        ),
        tx2gene="salmon/tx2gene.tsv",
    output:
        tsv="salmon/Raw.genes.tsv",
    message:
        "Aggregating genes counts with their gene names"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 4,
        time_min=lambda wildcards, attempt: attempt * 15,
        tmpdir="tmp",
    params:
        header=False,
        position=False,
        gencode=True,
        genes=True,
        index_label=True,
        fillna="Unknown",
        column="NumReads",
    log:
        "logs/aggregate/genes.log",
    wrapper:
        "bio/pandas/salmon"


# Aggregate per-gene TPM counts
# Usefull for GSEA, PCA, ...
rule aggregate_gene_counts:
    input:
        quant=expand(
            "salmon/pseudo_mapping/{sample}/quant.genes.sf", sample=design["Sample_id"]
        ),
        tx2gene="salmon/tx2gene.tsv",
    output:
        tsv="salmon/TPM.genes.tsv",
    message:
        "Aggregating genes counts with their gene names"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 4,
        time_min=lambda wildcards, attempt: attempt * 15,
        tmpdir="tmp",
    params:
        header=False,
        position=False,
        gencode=True,
        genes=True,
        index_label=True,
        fillna="Unknown",
    log:
        "logs/aggregate/genes.log",
    wrapper:
        "bio/pandas/salmon"


# Aggregate transcript counts
# Usefull for PCA, sample correlations, etc.
rule aggregate_transcript_counts:
    input:
        quant=expand(
            "salmon/pseudo_mapping/{sample}/quant.sf", sample=design["Sample_id"]
        ),
        tx2gene="salmon/tx2gene.tsv",
    output:
        tsv="salmon/TPM.transcripts.tsv",
    message:
        "Aggregating transcript counts with their gene names"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 4,
        time_min=lambda wildcards, attempt: attempt * 15,
        tmpdir="tmp",
    params:
        header=False,
        position=False,
        gencode=True,
        genes=False,
        index_label=True,
    log:
        "logs/aggregate/transcripts.log",
    wrapper:
        "bio/pandas/salmon"


# Required to annotate all previously described
# count tables with human-readable gene names
rule tx_to_gene:
    input:
        gtf=config["ref"]["gtf"],
    output:
        tx2gene=temp("salmon/tx2gene.tsv"),
        tx2gene_large=temp("salmon/tx2gene_with_positions.tsv"),
        gene2gene=temp("salmon/gene2gene.tsv"),
        gene2gene_large=temp("salmon/gene2gene_with_chr.tsv"),
    message:
        "Gathering transcripts and genes names together from GTF"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 4,
        time_min=lambda wildcards, attempt: min(attempt * 10, 15),
        tmpdir="tmp",
    log:
        "logs/tx_to_gene.log",
    wrapper:
        "bio/gtf/tx2gene"
