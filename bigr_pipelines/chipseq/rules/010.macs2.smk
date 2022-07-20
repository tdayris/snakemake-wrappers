rule macs2_callpeak_narrow:
    input:
        treatment="samtools/view/{sample}.bam",
        treatment_idx="samtools/view/{sample}.bam.bai",
    output:
        multiext(
            "macs2/callpeak/narrowPeak/{sample}",
            "_peaks.xls",
            "_peaks.narrowPeak",
            "_summits.bed",
        ),
    threads: 1
    resources:
        mem_mb=get_4gb_per_attempt,
        time_min=get_1h_per_attempt,
        tmpdir="tmp",
    log:
        "macs2/callpeak/narrowPeak/{sample}.log",
    params:
        "-g 2652783500 -f BAMPE",
    wrapper:
        "bio/macs2/callpeak"


rule macs2_callpeak_broad:
    input:
        treatment="samtools/view/{sample}.bam",
        treatment_idx="samtools/view/{sample}.bam.bai",
    output:
        multiext(
            "macs2/callpeak/broadPeak/{sample}",
            "_peaks.xls",
            "_treat_pileup.bdg",
            "_control_lambda.bdg",
            "_peaks.broadPeak",
            "_peaks.gappedPeak",
        ),
    threads: 1
    resources:
        mem_mb=get_4gb_per_attempt,
        time_min=get_1h_per_attempt,
        tmpdir="tmp",
    log:
        "logs/macs2/callpeak/broadPeak/{sample}.log",
    params:
        "-g 2652783500 -f BAMPE",
    wrapper:
        "bio/macs2/callpeak"


rule macs2_to_bed:
    input:
        "macs2/callpeak/{peaktype}/{sample}_peaks.{peaktype}",
    output:
        protected("macs2/callpeak/{peaktype}/{sample}_peaks.{peaktype}.bed"),
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir="tmp",
    log:
        "logs/macs2/tobed/{sample}.{peaktype}.log",
    params:
        "'BEGIN{FS=\"\\t\"} {print $1 FS $2 FS $3 FS $4 FS . FS +}'",
    conda:
        "../envs/bash.yaml"
    shell:
        "awk {params} {input} > {output} 2> {log}"
