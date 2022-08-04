rule multiqc_variant:
    input:
        html=expand(
            "fastp/html/pe/{sample}_{status}.fastp.html",
            sample=design["Sample_id"],
            status=status_list,
        ),
        json=expand(
            "fastp/json/pe/{sample}_{status}.fastp.json",
            sample=design["Sample_id"],
            status=status_list,
        ),
        sambamba_metrics=expand(
            "sambamba/markdup/{sample}_{status}.bam",
            sample=design["Sample_id"],
            status=status_list,
        ),
        fastq_screen=expand(
            "fastq_screen/{sample}.{stream}.{status}.fastq_screen.{ext}",
            sample=design["Sample_id"],
            stream=streams,
            ext=["txt", "png"],
            status=status_list,
        ),
        samtools_stats=expand(
            "samtools/stats/{sample}_{status}.{cleaning}.stats",
            sample=sample_list,
            status=status_list,
            cleaning=cleaning_status,
        ),
        csvstats=expand("snpeff/csvstats/{sample}.csv", sample=sample_list),
    output:
        protected("data_output/MultiQC/Somatic_Variant_Calling.html"),
    threads: 1
    resources:
        mem_mb=get_2gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir="tmp",
    retries: 1
    log:
        "logs/multiqc.log",
    wrapper:
        "bio/multiqc"


rule multiqc_mapping:
    input:
        html=expand(
            "fastp/html/pe/{sample}_{status}.fastp.html",
            sample=design["Sample_id"],
            status=status_list,
        ),
        json=expand(
            "fastp/json/pe/{sample}_{status}.fastp.json",
            sample=design["Sample_id"],
            status=status_list,
        ),
        sambamba_metrics=expand(
            "sambamba/markdup/{sample}_{status}.bam",
            sample=design["Sample_id"],
            status=status_list,
        ),
        fastq_screen=expand(
            "fastq_screen/{sample}.{stream}.{status}.fastq_screen.{ext}",
            sample=design["Sample_id"],
            stream=streams,
            ext=["txt", "png"],
            status=status_list,
        ),
        samtools_stats=expand(
            "samtools/stats/{sample}_{status}.{cleaning}.stats",
            sample=sample_list,
            status=status_list,
            cleaning=cleaning_status,
        ),
    output:
        protected("data_output/MultiQC/MappingQC.html"),
    threads: 1
    resources:
        mem_mb=get_2gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir="tmp",
    retries: 1
    log:
        "logs/multiqc.log",
    wrapper:
        "bio/multiqc"
