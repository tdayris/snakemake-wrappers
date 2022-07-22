index_datasets_config = {"genome": config["ref"]["fasta"]}


module index_datasets:
    snakefile:
        str(
            workflow_source_dir
            / ".."
            / ".."
            / "meta"
            / "bio"
            / "index_datasets"
            / "test"
            / "Snakefile"
        )
    config:
        index_datasets_config


use rule samtools_faidx from index_datasets


use rule picard_create_sequence_dictionnary from index_datasets


rule sambamba_index_bam:
    input:
        "{tool}/{subcommand}/{sample}_{status}.bam",
    output:
        "{tool}/{subcommand}/{sample}_{status}.bam.bai",
    message:
        "Indexing {wildcards.sample} ({wildcards.status}) "
        "from {wildcards.tool}:{wildcards.subcommand}"
    threads: 8
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 16,
        time_min=lambda wildcards, attempt: attempt * 35,
        tmpdir="tmp",
    log:
        "sambamba/index/{tool}_{subcommand}/{sample}_{status}.log",
    params:
        extra="",
    wrapper:
        "bio/sambamba/index"
