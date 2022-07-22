rule deeptools_bamcov:
    input:
        bam="samtools/view/{sample}.bam",
        bai="samtools/view/{sample}.bam.bai",
        blacklist=config["reference"]["blacklist"],
    output:
        protected("deeptools/bamcoverage/{sample}.bw"),
    threads: 20
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_1h_per_attempt,
        tmpdir="tmp",
    log:
        "logs/deeptools/bamcov/{sample}.log",
    conda:
        "../envs/deeptools.yaml"
    params:
        extra=(
            "--binSize 10 "
            "--effectiveGenomeSize 2652783500 "
            "--normalizeUsing RPKM "
            "--skipNonCoveredRegions "
            "--extendReads "
            "--centerReads "
        ),
    shell:
        "bamCoverage "
        "--bam {input.bam} "
        "--outFileName {output} "
        "--outFileFormat 'bigwig' "
        "--blackListFileName {input.blacklist} "
        "{params.extra} "
        "> {log} 2>&1"


rule deeptools_fingerprint:
    input:
        bam_files=expand("samtools/view/{sample}.bam", sample=sample_list),
        bam_idx=expand("samtools/view/{sample}.bam.bai", sample=sample_list),
        blacklist=config["reference"]["blacklist"],
    output:
        fingerprint="deeptools/plot_fingerprint/plot_fingerprint.png",
        counts=temp("deeptools/plot_fingerprint/raw_counts.tab"),
        qc_metrics=temp("deeptools/plot_fingerprint/qc_metrics.txt"),
    threads: 20
    resources:
        mem_mb=get_2gb_per_attempt,
        time_min=get_1h_per_attempt,
        tmpdir="tmp",
    params:
        "--extendReads --centerReads --skipZeros --verbose ",
    wrapper:
        "bio/deeptools/plotfingerprint"


rule deeptools_compute_matrix:
    input:
        bed=expand(
            "macs2/callpeak/{peaktype}/{sample}_peaks.{peaktype}.bed",
            sample=sample_list,
            allow_missing=True,
        ),
        bigwig=expand(
            "deeptools/bamcoverage/{sample}.bw",
            sample=sample_list,
        ),
        blacklist=config["reference"]["blacklist"],
    output:
        matrix_gz=temp("deeptools/matrix_files/{peaktype}.gz"),
    threads: 20
    resources:
        mem_mb=get_4gb_per_attempt,
        time_min=get_2h_per_attempt,
        tmpdir="tmp",
    log:
        "logs/deeptools/computematrix/{peaktype}.log",
    params:
        extra=lambda wildcards, input: "--verbose --numberOfProcessors 20 --blackListFileName {input.blacklist} --skipZeros",
        command="reference-point",
    wrapper:
        "bio/deeptools/computematrix"


rule deeptools_plot_heatmap:
    input:
        "deeptools/matrix_files/{peaktype}.gz",
    output:
        heatmap_img="deeptools/plot_heatmap/{peaktype}.heatmap.png",
        regions=temp("deeptools/plot_heatmap/{peaktype}.heatmap_regions.bed"),
        heatmap_matrix=temp("deeptools/plot_heatmap/{peaktype}.heatmap_matrix.tab"),
    threads: 20
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_1h_per_attempt,
        tmpdir="tmp",
    log:
        "logs/deeptools/heatmap/{peaktype}.log",
    params:
        " --kmeans ",
    wrapper:
        "bio/deeptools/plotheatmap"


rule deeptools_multibigwigsummary:
    input:
        bins=expand("deeptools/bamcoverage/{sample}.bw", sample=sample_list),
        blacklist=config["reference"]["blacklist"],
    output:
        temp("deeptools/multibigwigsummary/results.npz"),
    threads: 20
    resources:
        mem_mb=get_4gb_per_attempt,
        time_min=get_2h_per_attempt,
        tmpdir="tmp",
    log:
        "logs/deeptools/multibigwigsummary.log",
    conda:
        "../envs/deeptools.yaml"
    params:
        "--verbose ",
    shell:
        "multiBigwigSummary bins "
        "--bwfiles {input.bins} "
        "--outFileName {output} "
        "--blackListFileName {input.blacklist} "
        "--numberOfProcessors {threads} "
        "{params} "
        "> {log} 2>&1 "


rule deeptools_plotcorrelation:
    input:
        cordata="deeptools/multibigwigsummary/results.npz",
    output:
        png="deeptools/correlations/sample_correlation.png",
        counts=temp("deeptools/correlations/SpearmanCorr_readCounts.tab"),
    threads: 1
    resources:
        mem_mb=get_4gb_per_attempt,
        time_min=get_1h_per_attempt,
        tmpdir="tmp",
    log:
        "logs/deeptools/correlation.log",
    conda:
        "../envs/deeptools.yaml"
    params:
        extra=(
            "--corMethod spearman "
            "--whatToPlot heatmap "
            "--skipZeros "
            "--plotFileFormat png "
        ),
    shell:
        "plotCorrelation "
        "--corData {input.cordata} "
        "--plotFile {output.png} "
        "--outFileCorMatrix {output.counts} "
        "{params.extra} "
        "> {log} 2>&1 "


rule deeptools_plotpca:
    input:
        cordata="deeptools/multibigwigsummary/results.npz",
    output:
        png="deeptools/pca/PCA.png",
        data="deeptools/pca/PCA.txt",
    threads: 1
    resources:
        mem_mb=get_4gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir="tmp",
    log:
        "logs/deeptools/pca.log",
    conda:
        "../envs/deeptools.yaml"
    params:
        extra="--plotFileFormat png",
    shell:
        "plotPCA "
        "--corData {input.cordata} "
        "--plotFile {output.png} "
        "--outFileNameData {output.data} "
        "> {log} 2>&1 "
