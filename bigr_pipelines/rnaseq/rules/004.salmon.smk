# This rule pseudo-map and quantifies your paired reads over the indexed
# reference.
"""
004.salmon_quant
from
-> 002.fastp_clean
by
-> 007.tximport
-> 005.subset_gene_counts
"""


rule salmon_quant:
    input:
        r1="002.fastp/trimmed/{sample}.1.fastq",
        r2="002.fastp/trimmed/{sample}.2.fastq",
        index=ancient(config["salmon"]["index"]),
        gtf=ancient(config["reference"]["gtf"]),
    output:
        quant=temp("004.salmon/pseudo_mapping/{sample}/quant.sf"),
        quant_genes=temp("004.salmon/pseudo_mapping/{sample}/quant.genes.sf"),
        lib=temp("004.salmon/pseudo_mapping/{sample}/lib_format_counts.json"),
        aux_info=temp(directory("004.salmon/pseudo_mapping/{sample}/aux_info")),
        cmd_info=temp("004.salmon/pseudo_mapping/{sample}/cmd_info.json"),
        libparams=temp(directory("004.salmon/pseudo_mapping/{sample}/libParams")),
        logs=temp(directory("004.salmon/pseudo_mapping/{sample}/logs")),
    threads: config.get("max_threads", 20)
    resources:
        time_min=get_salmon_time_per_attempt,
        mem_mb=get_15gb_per_attempt,
        tmpdir="tmp",
    retries: 1
    params:
        libType=config["salmon"].get("salmon_libtype", "A"),
        extra=config["salmon"].get(
            "salmon_quant_extra", "--numBootstraps 100 --gcBias --seqBias --posBias"
        ),
    log:
        "logs/004.salmon/quant/{sample}.log",
    wrapper:
        "bio/salmon/quant"


# Required to annotate all previously described
# count tables with human-readable gene names
"""
004.tx_to_gene
from
-> Entry job
by
-> 004.aggregate_raw_counts
-> 004.aggregate_gene_counts
-> 007.tximport
"""


rule tx_to_gene:
    input:
        gtf=config["reference"]["gtf"],
    output:
        tx2gene_small=temp("004.salmon/tx2gene_small.tsv"),
        tx2gene=temp("004.salmon/tx2gene.tsv"),
        tx2gene_large=temp("004.salmon/tx2gene_with_positions.tsv"),
        gene2gene=temp("004.salmon/gene2gene.tsv"),
        gene2gene_large=temp("004.salmon/gene2gene_with_chr.tsv"),
    threads: 1
    resources:
        mem_mb=get_4gb_per_gb_input_per_attempt,
        time_min=get_15min_per_gb_input_per_attempt,
        tmpdir="tmp",
    retries: 1
    log:
        "logs/004.tx_to_gene.log",
    wrapper:
        "bio/gtf/tx2gene"


# Work on raw counts, usefull for manual
# DESeq2/EdgeR/...
"""
004.aggregate_raw_counts
from
-> 004.tx_to_gene
-> 004.salmon_quant
by
-> End job
"""


rule aggregate_raw_counts:
    input:
        quant=expand(
            "004.salmon/pseudo_mapping/{sample}/quant.genes.sf",
            sample=design["Sample_id"],
        ),
        tx2gene="004.salmon/tx2gene.tsv",
    output:
        tsv=protected("data_output/Quantification/Raw.genes.tsv"),
    threads: 1
    resources:
        mem_mb=get_4gb_per_gb_input_per_attempt,
        time_min=get_15min_per_gb_input_per_attempt,
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
        "logs/004.aggregate/genes.log",
    wrapper:
        "bio/pandas/salmon"


# Aggregate per-gene TPM counts
# Usefull for GSEA, PCA, ...
"""
004.aggregate_gene_counts
from
-> 004.tx_to_gene
-> 004.salmon_quant
by
-> 005.subset_gene_counts
"""


rule aggregate_gene_counts:
    input:
        quant=expand(
            "004.salmon/pseudo_mapping/{sample}/quant.genes.sf",
            sample=design["Sample_id"],
        ),
        tx2gene="004.salmon/tx2gene.tsv",
    output:
        tsv=protected("data_output/Quantification/TPM.genes.tsv"),
    threads: 1
    resources:
        mem_mb=get_4gb_per_gb_input_per_attempt,
        time_min=get_15min_per_gb_input_per_attempt,
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
        "logs/004.aggregate/genes.log",
    wrapper:
        "bio/pandas/salmon"


# Aggregate transcript counts
# Usefull for PCA, sample correlations, etc.
"""
004.aggregate_transcript_counts
from
-> 004.tx_to_gene
-> 004.salmon_quant
by
-> End job
"""


rule aggregate_transcript_counts:
    input:
        quant=expand(
            "004.salmon/pseudo_mapping/{sample}/quant.sf", sample=design["Sample_id"]
        ),
        tx2gene="004.salmon/tx2gene.tsv",
    output:
        tsv=protected("data_output/Quantification/TPM.transcripts.tsv"),
    threads: 1
    resources:
        mem_mb=get_4gb_per_gb_input_per_attempt,
        time_min=get_15min_per_gb_input_per_attempt,
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
        "logs/004.aggregate/transcripts.log",
    wrapper:
        "bio/pandas/salmon"
