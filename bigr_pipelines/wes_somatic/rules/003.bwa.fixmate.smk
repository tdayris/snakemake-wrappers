"""
This rule sorts reads by position for further analyses.
This rule is shadowed in order to automatically delete undetermined number of
temporary files on error, which leads retry to fail
"""


rule sambamba_sort_coordinate:
    input:
        mapping="samtools/fixmate/{sample}_{status}.bam",
    output:
        mapping=temp("sambamba/sort/{sample}_{status}.bam"),
    threads: min(config.get("max_threads", 20), 10)
    resources:
        mem_mb=get_4gb_per_attempt,
        time_min=get_90min_per_attempt,
        tmpdir=tmp,
    retries: 1
    shadow:
        "shallow"
    params:
        extra=config["sambamba"].get("sort_extra", ""),
    log:
        "logs/sambamba/sort/{sample}.{status}.log",
    wrapper:
        "bio/sambamba/sort"


rule sambamba_sort_raw:
    input:
        mapping="bwa_mem2/mem/{sample}_{status}.bam",
    output:
        mapping=temp("bwa_mem2/sorted/{sample}_{status}.bam"),
    threads: min(config.get("max_threads", 20), 10)
    resources:
        mem_mb=get_2gb_per_attempt,
        time_min=get_90min_per_attempt,
        tmpdir=tmp,
    retries: 1
    shadow:
        "shallow"
    params:
        extra=config["sambamba"].get("sort_extra", ""),
    log:
        "logs/bwa/sort/{sample}.{status}.log",
    wrapper:
        "bio/sambamba/sort"


"""
BWA sometimes fails to annotate read mates correctly. We fix this behaviour
with the rule below.
"""


rule samtools_fixmate:
    input:
        "bwa_mem2/mem/{sample}_{status}.bam",
    output:
        temp("samtools/fixmate/{sample}_{status}.bam"),
    threads: config.get("max_threads", 20)
    resources:
        mem_mb=get_4gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir=tmp,
    retries: 1
    params:
        config["samtools"].get("fixmate_extra", "-cmr"),
    log:
        "logs/samtools/fixmate/{sample}.{status}.log",
    wrapper:
        "bio/samtools/fixmate"


rule samtools_view_bwa:
    input:
        "bwa_mem2/mem/{sample}_{status}.sam",
        fasta=config["reference"]["fasta"],
        fasta_idx=get_fai(config["reference"]["fasta"]),
        fasta_dict=get_dict(config["reference"]["fasta"]),
        bed=config["reference"]["capture_kit_bed"],
    output:
        temp("bwa_mem2/mem/{sample}_{status}.unsorted.bam"),
    threads: min(config.get("max_threads", 20), 10)
    resources:
        mem_mb=get_2gb_per_attempt,
        time_min=get_35min_per_attempt,
        tmpdir=tmp,
    retries: 1
    params:
        extra="-h",
    log:
        "logs/sambamba/view/{sample}.{status}.raw_star.log",
    wrapper:
        "bio/samtools/view"


rule sambamba_sort_raw_bwa:
    input:
        mapping="bwa_mem2/mem/{sample}_{status}.unsorted.bam",
    output:
        mapping=temp("bwa_mem2/mem/{sample}_{status}.bam"),
    threads: min(config.get("max_threads", 20), 10)
    resources:
        mem_mb=get_2gb_per_attempt,
        time_min=get_90min_per_attempt,
        tmpdir=tmp,
    retries: 1
    shadow:
        "shallow"
    params:
        extra="--sort-by-name",
    log:
        "logs/star/sort/{sample}.{status}.log",
    wrapper:
        "bio/sambamba/sort"


"""
This rule maps your reads against the indexed reference with BWA.
"""


rule bwa_mem:
    input:
        reads=expand(
            "fastp/trimmed/{sample}_{status}.{stream}.fastq",
            stream=["1", "2"],
            allow_missing=True,
        ),
        idx=config["reference"]["bwa_index"],
    output:
        temp("bwa_mem2/mem/{sample}_{status}.sam"),
    threads: config.get("max_threads", 20)
    resources:
        mem_mb=get_75gb_and_2gb_per_attempt,
        time_min=get_2h_per_attempt,
        tmpdir=tmp,
    retries: 1
    shadow:
        "shallow"
    params:
        index=lambda wildcards, input: os.path.splitext(input["index"][0])[0],
        extra=r"-R '@RG\tID:{sample}_{status}\tSM:{sample}_{status}\tPU:{sample}_{status}\tPL:ILLUMINA\tCN:IGR\tDS:WES\tPG:BWA-MEM2' -M -A 2 -E 1",
        sort="samtools",
        sort_order="queryname",
        sort_extra="-m 1536M",
    log:
        "logs/bwa_mem2/mem/{sample}.{status}.log",
    wrapper:
        "bio/bwa-mem2/mem"
