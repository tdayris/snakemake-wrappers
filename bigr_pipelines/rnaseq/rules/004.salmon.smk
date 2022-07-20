# This rule pseudo-map and quantifies your paired reads over the indexed
# reference.
rule salmon_quant_paired:
    input:
        r1="fastp/trimmed/{sample}.1.fastq",
        r2="fastp/trimmed/{sample}.2.fastq",
        index=ancient(config["salmon"]["index"]),
        gtf=ancient(config["reference"]["gtf"]),
    output:
        quant="salmon/pseudo_mapping/{sample}/quant.sf",
        lib="salmon/pseudo_mapping/{sample}/lib_format_counts.json",
    threads: config.get("max_threads", 20)
    resources:
        time_min=get_45min_per_attempt,
        mem_mb=get_10gb_per_attempt,
        tmpdir="tmp",
    retries: 2
    params:
        libType=config["salmon"].get("salmon_libtype", "A"),
        extra=config["salmon"].get(
            "salmon_quant_extra", "--numBootstraps 100 --gcBias --seqBias --posBias"
        ),
    log:
        "logs/salmon/quant/{sample}.log",
    wrapper:
        "bio/salmon/quant"


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
    threads: 1
    resources:
        mem_mb=get_4gb_per_attempt,
        time_min=get_15min_per_attempt,
        tmpdir="tmp",
    retries: 1
    log:
        "logs/tx_to_gene.log",
    wrapper:
        "bio/gtf/tx2gene"


# Work on raw counts, usefull for manual
# DESeq2/EdgeR/...
rule aggregate_raw_counts:
    input:
        quant=expand(
            "salmon/pseudo_mapping/{sample}/quant.genes.sf", sample=design["Sample_id"]
        ),
        tx2gene="salmon/tx2gene.tsv",
    output:
        tsv=report(
            "data_output/Quantification/Raw.genes.tsv",
            caption=str(worflow_source_dir / "reports" / "salmon.raw.genes.rst"),
            category="Counts",
        ),
    threads: 1
    resources:
        mem_mb=get_4gb_per_attempt,
        time_min=get_15min_per_attempt,
        tmpdir="tmp",
    retries: 1
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
        tsv=report(
            "data_output/Quantification/TPM.genes.tsv",
            caption=str(worflow_source_dir / "reports" / "salmon.tpm.genes.rst"),
            category="Counts",
        ),
    threads: 1
    resources:
        mem_mb=get_4gb_per_attempt,
        time_min=get_15min_per_attempt,
        tmpdir="tmp",
    retries: 1
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
        tsv=report(
            "data_output/Quantification/TPM.transcripts.tsv",
            caption=str(worflow_source_dir / "reports" / "salmon.tpm.transcripts.rst"),
            category="Counts",
        ),
    threads: 1
    resources:
        mem_mb=get_4gb_per_attempt,
        time_min=get_15min_per_attempt,
        tmpdir="tmp",
    retries: 1
    params:
        header=False,
        position=False,
        gencode=True,
        genes=False,
        index_label=True,
        fillna="Unknown",
    log:
        "logs/aggregate/transcripts.log",
    wrapper:
        "bio/pandas/salmon"
