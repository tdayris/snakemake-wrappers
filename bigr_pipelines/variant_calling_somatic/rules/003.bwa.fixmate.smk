module bwa_fixmate_meta:
    snakefile:
        str(
            workflow_source_dir
            / ".."
            / ".."
            / "meta"
            / "bio"
            / "bwa_fixmate"
            / "test"
            / "Snakefile"
        )
    config:
        {"threads": config["threads"], "genome": config["ref"]["fasta"]}


use rule sambamba_index from bwa_fixmate_meta with:
    input:
        "sambamba/sort/{sample}_{status}.bam",
    output:
        temp("sambamba/sort/{sample}_{status}.bam.bai"),
    message:
        "Indexing mapped reads of {wildcards.status} {wildcards.sample}"
    log:
        "logs/sambamba/sort/{sample}.{status}.log",


use rule sambamba_sort_coordinate from bwa_fixmate_meta with:
    input:
        mapping="samtools/fixmate/{sample}_{status}.bam",
    output:
        mapping=temp("sambamba/sort/{sample}_{status}.bam"),
    message:
        "Sorting {wildcards.status} {wildcards.sample} reads by position"
    log:
        "logs/sambamba/sort/{sample}.{status}.log",


use rule samtools_fixmate from bwa_fixmate_meta with:
    input:
        "bwa_mem2/mem/{sample}_{status}.bam",
    output:
        temp("samtools/fixmate/{sample}_{status}.bam"),
    message:
        "Fixing mate annotation on {wildcards.status} "
        "{wildcards.sample} with Samtools"
    log:
        "logs/samtools/fixmate/{sample}.{status}.log",


use rule bwa_mem from bwa_fixmate_meta with:
    input:
        reads=expand(
            "fastp/trimmed/pe/{sample}_{status}.{stream}.fastq",
            stream=["1", "2"],
            allow_missing=True,
        ),
        index=config["ref"]["bwa_index"],
    output:
        temp("bwa_mem2/mem/{sample}_{status}.bam"),
    message:
        "Mapping {wildcards.status} {wildcards.sample} with BWA"
    params:
        index=lambda wildcards, input: os.path.splitext(input["index"][0])[0],
        extra="-R '@RG\tID:{sample}_{status}\tSM:{sample}_{status}\tPU:{sample}_{status}\tPL:ILLUMINA\tCN:IGR\tDS:WES\tPG:BWA-MEM2' -M -A 2 -E 1",
        sort="samtools",  # We chose Samtools to sort by queryname
        sort_order="queryname",  # Queryname sort is needed for a fixmate
        sort_extra="-m 1536M",  # We extand the sort buffer memory
    log:
        "logs/bwa_mem2/mem/{sample}.{status}.log",
