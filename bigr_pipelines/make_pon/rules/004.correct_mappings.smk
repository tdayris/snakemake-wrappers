"""
004::samtools_fixmate:
-> 003::bwa_mem_align

BWA sometimes fails to annotate read mates correctly. We fix this behaviour
with the rule below.
"""
rule samtools_fixmate:
    input:
        "bwa/mem/{sample}.bam"
    output:
        temp("samtools/fixmate/{sample}.bam")
    threads: min(config.get("max_threads", 20), 20)
    resources:
        mem_mb=get_4gb_per_attempt,
        time_min=get_75min_per_attempt,
        tmpdir="tmp"
    params:
        config["params"].get("fixmate_extra", "-c -m -r")
    log:
        "logs/samtools/fixmate/{sample}.log"
    wrapper:
        "bio/samtools/fixmate"


"""
004::sambamba_sort_coordinate:
-> 004::samtools_fixmate

This rule sorts reads by position for further analyses.
This rule is shadowed in order to automatically delete undetermined number of
temporary files on error, which leads retry to fail
"""
rule sambamba_sort_coordinate:
    input:
        mapping="samtools/fixmate/{sample}.bam"
    output:
        mapping=temp("sambamba/sort/{sample}.bam"),
    threads: min(config.get("max_threads", 20), 8)
    resources:
        mem_mb=get_4gb_per_attempt,  # Make sure memory matches --memory-limit
        time_min=get_1h_per_attempt,
        tmpdir="tmp"
    shadow: "shallow"
    params:
        extra = "--show-progress"
    log:
        "logs/sambamba/sort/{sample}.log"
    wrapper:
        "bio/sambamba/sort"


"""
004::sambamba_markdup
-> 004::sambamba_sort_coordinate

This rule removes deuplicate reads
"""
rule sambamba_markdup:
    input:
        "sambamba/sort/{sample}.bam"
    output:
        "sambamba/markdup/{sample}.bam"
    resources:
        mem_mb=get_4gb_per_attempt,
        time_min=get_1h_per_attempt,
        tmpdir="tmp"
    log:
        "logs/sambamba/markdup.log"
    params:
        extra=config["params"].get("sambamba_markdup", "--remove-duplicates --overflow-list-size 600000")
    wrapper:
        f"{wrappers}/bio/sambamba/markdup"


"""
004::samtools_view_filter:
-> 004::sambamba_markdup

This rule filters the BAM file over capture kit bed
"""
rule samtools_view_filter:
    input:
        bam="sambamba/markdup/{sample}.bam",
        bed=config[genome_id]["bed"],
        ref=config[genome_id]["fasta"],
        ref_idx=fai_file
    output:
        bam=temp("samtools/filter/{sample}.bam")
    threads: 5
    resources:
        mem_mb=get_1p5gb_per_attempt,
        time_min=get_20min_per_attempt,
        tmpdir="tmp"
    log:
        "logs/samtools/view/{sample}.filter.log"
    params:
        "-h -b"
    conda:
        "envs/samtools.yaml"
    shell:
        "samtools view {params} " # Tool and optional parameters
        "-T {input.ref} " # Path to reference genome sequence
        "-t {input.ref_idx} " # Path to reference genome index
        "--region-file {input.bed} " # Path to capture kit bed
        "-o {output.bam} " # Path to output cram file
        "{input.bam} " # Path to input bam file
        "> {log} 2>&1" 


"""
004::sambamba_index:
-> 004::samtools_view_filter

This rule indexes the bam file since almost all downstream tools requires it
"""
rule sambamba_index:
    input:
        "samtools/filter/{sample}.bam"
    output:
        temp("samtools/filter/{sample}.bam.bai")
    threads: min(config.get("max_threads", 20), 5)
    resources:
        mem_mb=get_1p5gb_per_attempt,
        time_min=get_20min_per_attempt,
        tmpdir="tmp"
    params:
        extra = "--show-progress"
    log:
        "logs/sambamba/index/{sample}.log"
    wrapper:
        "bio/sambamba/index"


"""
004::samtools_view_cram:
-> 004::samtools_view_filter
-> 004::sambamba_index
-> 001::samtools_index_genome

This rule compress BAM files in CRAM for future storage
"""
rule samtools_view_cram:
    input:
        bam="samtools/filter/{sample}.bam",
        bam_idx="samtools/filter/{sample}.bam.bai",
        ref=config[genome_id]["fasta"],
        ref_idx=ancient(fai_file),
    output:
        cram="mapping/{sample}.cram",
    resources:
        mem_mb=get_1p5gb_per_attempt,
        time_min=get_1h_per_attempt,
        tmpdir="tmp"
    log:
        "logs/samtools/view/{sample}.cram.log"
    params:
        "-h -C"
    conda:
        "envs/samtools.yaml"
    shell:
        "samtools view {params} " # Tool and optional parameters
        "-T {input.ref} " # Path to reference genome sequence
        "-t {input.ref_idx} " # Path to reference genome index
        "-o {output.cram} " # Path to output cram file
        "{input.bam} " # Path to input bam file
        "> {log} 2>&1"
