"""
007::multiqc_report:
-> 002::fastp_trimming
-> 007::picard_collectmultiplemetrics

Aggregates quality control reports
"""
rule multiqc_report:
    input:
        fastp_json=expand(
            "fastp/json/pe/{sample}.fastp.json",
            sample=design.index
        ),
        fastp_html=expand(
            "fastp/html/pe/{sample}.fastp.html",
            sample=design.index
        ),
        mappings=expand(
            multiext(
                "picard/collectmultiplemetrics/{sample}/{sample}",
                ".alignment_summary_metrics",
                ".insert_size_metrics",
                ".insert_size_histogram.pdf",
                ".gc_bias.detail_metrics",
                ".gc_bias.summary_metrics",
                ".gc_bias.pdf",
                ".bait_bias_detail_metrics",
                ".bait_bias_summary_metrics",
                ".error_summary_metrics",
                ".pre_adapter_detail_metrics",
                ".pre_adapter_summary_metrics",
                ".quality_distribution_metrics",
                ".quality_distribution.pdf",
            ),
            sample=design.index
        ),
        bcftools=expand(
            "bcftools/stats/{sample}.stats.txt",
            sample=design.index
        )
    output:
        "../results/multiqc/PoN.html"
    threads: 1
    resources:
        mem_mb=get_1p5gb_per_attempt,
        time_min=get_10_minutes_per_gb,
        tmpdir="tmp",
    log:
        "logs/sambamba/bwa/{sample}.log"
    params:
        extra=""
    wrapper:
        "bio/multiqc"



"""
007::sambamba_sort_coordinate_raw_bam:
-> 003::bwa_mem_align

Sort raw bwa bam file in order to gather QC
"""
rule sambamba_sort_coordinate_raw_bam:
    input:
        "../results/bwa/mem/{sample}.bam"
    output:
        temp("../results/bwa/coordinates/{sample}.bam")
    threads: min(config.get("max_threads", 20), 20)
    resources:
        mem_mb=get_1p5gb_per_attempt,
        time_min=get_10_minutes_per_gb,
        tmpdir="tmp",
    log:
        "logs/sambamba/bwa/{sample}.log"
    params:
        extra=""
    wrapper:
        "bio/sambamba/sort"


"""
007::sambamba_index_coordinate_raw_bam
-> 007::sambamba_sort_coordinate_raw_bam

Index sorted raw bam
"""
rule sambamba_index_coordinate_raw_bam:
    input:
        "../results/bwa/coordinates/{sample}.bam"
    output:
        temp("../results/bwa/coordinates/{sample}.bam.bai")
    resources:
        mem_mb=get_1p5gb_per_attempt,
        time_min=get_10_minutes_per_attempt,
        tmpdir="tmp",
    log:
        "logs/sambamba/bwa/{sample}.index.log"
    params:
        extra=""
    wrapper:
        "bio/sambamba/index"
    

"""
007::picard_collectmultiplemetrics:
-> 007::sambamba_index_coordinate_raw_bam
-> 007::sambamba_sort_coordinate_raw_bam
-> 001::samtools_index_genome
-> 001::picard_createsequencedictionary
"""
rule picard_collect_multiple_metrics:
    input:
        bam="../results/bwa/coordinates/{sample}.bam",
        bai="../results/bwa/coordinates/{sample}.bam.bai",
        ref=fasta_file,
        ref_idx=fai_file,
        ref_dict=dict_file,
    output:
        multiext(
            "picard/collectmultiplemetrics/{sample}/{sample}",
            ".alignment_summary_metrics",
            ".insert_size_metrics",
            ".insert_size_histogram.pdf",
            ".gc_bias.detail_metrics",
            ".gc_bias.summary_metrics",
            ".gc_bias.pdf",
            ".bait_bias_detail_metrics",
            ".bait_bias_summary_metrics",
            ".error_summary_metrics",
            ".pre_adapter_detail_metrics",
            ".pre_adapter_summary_metrics",
            ".quality_distribution_metrics",
            ".quality_distribution.pdf",
        ),
    resources:
        mem_mb=get_4gb_per_gb,
        time_min=get_1h_per_attempt,
        tmpdir="tmp"
    log:
        "logs/picard/collectmultiplemetrics/{sample}.log"
    params:
        extra=(
            "--VALIDATION_STRINGENCY LENIENT "
            "--METRIC_ACCUMULATION_LEVEL null "
            "--METRIC_ACCUMULATION_LEVEL SAMPLE"
        )
    wrapper:
        "bio/picard/collectmultiplemetrics"


rule bcftools_stats:
    input:
        ref=config[genome_id]["fasta"],
        ref_idx=fai_file,
        region=config[genome_id]["bed"],
        vcf="mutect2/filter/{sample}.vcf.gz"
        vcf_idx="mutect2/filter/{sample}.vcf.gz.tbi"
    output:
        temp("bcftools/stats/{sample}.stats.txt")
    threads: 1
    resources:
        mem_mb=get_1p5gb_per_attempt,
        time_min=get_1h_per_attempt,
        tmpdir="tmp"
    conda:
        str(workflow_source_dir / "envs" / "bcftools.yaml")
    params:
        extra=""
    script:
        str(workflow_source_dir / "scripts" / "007.bcftools_stats.py")
