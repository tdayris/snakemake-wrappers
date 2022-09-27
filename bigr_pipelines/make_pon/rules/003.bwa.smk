"""
002::bwa_mem_index Optional rule
-> 001::samtools_index_genome

Index genome sequence, in case user did not
provided any gneome index in configfile
"""


rule bwa_mem_index:
    input:
        config[genome_id]["fasta"],
        ancient(fai_file),
    output:
        bwa_sequence_index,
    threads: min(config.get("max_threads", 20), 20)
    resources:
        mem_mb=get_75gb_and_2gb_per_attempt,
        time_min=get_4h_per_attempt,
        tmpdir="tmp",
    log:
        "logs/bwa/index.log",
    wrapper:
        "bio/bwa-mem2/index"


"""
003::bwa_mem_align:
-> 002::fastp_trimming
-> 003::bwa_mem_index

This rule aligns reads to the genome sequences

This rule should always be used in the pipeline
"""


rule bwa_mem_align:
    input:
        idx=ancient(bwa_sequence_index),
        reads=expand(
            "fastp/trimmed/{sample}.{stream}.fastq",
            stream=["1", "2"],
            allow_missing=True,
        ),
    output:
        temp("bwa/mem/{sample}.bam"),
    threads: min(config.get("max_threads", 20), 20)
    resources:
        mem_mb=get_75gb_and_2gb_per_attempt,
        time_min=get_4h_per_attempt,
        tmpdir="tmp",
    shadow:
        "shallow"  # When samtools sort fails, tmp files are not deleted
    log:
        "logs/bwa/mem/{sample}.log",
    params:
        index=lambda wildcards, input: os.path.splitext(input["index"][0])[0],
        extra=r"-R '@RG\tID:{sample}_{status}\tSM:{sample}_{status}\tPU:{sample}_{status}\tPL:ILLUMINA\tCN:IGR\tDS:WES\tPG:BWA-MEM2' -M -A 2 -E 1",
        sort="samtools",
        sort_order="queryname",
        sort_extra="-m 1536M",
    wrapper:
        "bio/bwa-mem2/mem"
