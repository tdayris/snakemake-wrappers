rule deeptools_bamcov:
    input:
        bam="samtools/view/{sample}.bam",
        bai="samtools/view/{sample}.bam.bai",
        blacklist=config["reference"]["blacklist"],
    output:
        protected("deeptools/bamcoverage/{sample}.bw")
    threads: 20
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_1h_per_attempt,
        tmpdir="tmp"
    log:
        "logs/deeptools/bamcov/{sample}.log"
    conda:
        "envs/deeptools.yaml"
    params:
        extra=(
            "--binSize 1 "
            "--effectiveGenomeSize 2652783500 "
            "--normalizeUsing RPKM "
            "--skipNonCoveredRegions "
            "--extendReads "
            "--centerReads "
        )
    shell:
        "bamCoverage "
        "--bam {input.bam} "
        "--outFileName {output} "
        "--outFileFormat 'bigwig' "
        "--blackListFileName {input.blacklist} "
        "{params.extra} "
        "> {log} 2>&1"


rule deeptools_compute_matrix:
    input:
        bed="macs2/callpeak/{peaktype}/{sample}_peaks.{peaktype}.bed",
        bigwig="deeptools/bamcoverage/{sample}.bw",
        blacklist=config["reference"]["blacklist"]
    output:
        matrix_gz=temp("deeptools/matrix_files/{sample}.gz"),
    threads: 20
    resources:
        mem_mb=get_4gb_per_attempt,
        time_min=get_2h_per_attempt,
        tmpdir="tmp",
    log:
        "logs/deeptools/computematrix/{sample}.log"
    params:
        extra=lambda wildcards, input:"--verbose --numberOfProcessors 20 --blackListFileName {input.blacklist} --skipZeros",
        command="reference-point",
    wrapper:
        "bio/deeptools/computematrix"


rule deeptools_plot_heatmap:
    input:
        "deeptools/matrix_files/{sample}.gz"
    output:
        heatmap_img="deeptools/plot_heatmap/heatmap.png",
    threads: 20
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_1h_per_attempt,
        tmpdir="tmp"
    log:
        "logs/deeptools/heatmap/{sample}.log"
    params:
        " --kmeans "
    wrapper:
        "bio/deeptools/plotheatmap"