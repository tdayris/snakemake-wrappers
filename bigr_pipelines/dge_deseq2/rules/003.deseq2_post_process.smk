deseq2_post_process_config = {
    "condition_dict": condition_dict,
    "samples_per_prefixes": samples_per_prefixes,
    "design": design.copy(),
    "thresholds": config["thresholds"],
    "genes_of_interest": config.get(
        "genes_of_interest", ["ENSG00000141510"]
    ),
    "chromosomes": config["ref"].get(
        "chromosomes",
        list(range(24)) + ["MT", "X", "Y"] + list(map(str, range(24)))
    )
}

snakefile_deseq2_post_process = snakemake.workflow.srcdir(
    "../../meta/bio/deseq2_post_process/test/Snakefile"
)


module deseq2_post_process:
    snakefile: snakefile_deseq2_post_process
    config: deseq2_post_process_config


use rule * from deseq2_post_process


use rule pandas_merge_salmon_tr from deseq2_post_process with:
    input:
        quant = expand(
            "salmon/pseudo_mapping/{sample}/quant.sf",
            sample=design.Sample_id.tolist()
        ),
        tx2gene = "tximport/transcripts2genes.tsv"


use rule multiqc from deseq2_post_process with:
    input:
        txt=lambda wildcards: expand(
            "fastq_screen/{sample}.{stream}.fastq_screen.txt",
            sample=samples_per_prefixes[wildcards.comparison],
            stream=["1", "2"]
        ),
        png=lambda wildcards: expand(
            "fastq_screen/{sample}.{stream}.fastq_screen.png",
            sample=samples_per_prefixes[wildcards.comparison],
            stream=["1", "2"]
        ),
        salmon=lambda wildcards: expand(
            "salmon/pseudo_mapping/{sample}/quant.sf",
            sample=samples_per_prefixes[wildcards.comparison]
        ),
        html=lambda wildcards: expand(
            "fastp/html/pe/{sample}.fastp.html",
            sample=samples_per_prefixes[wildcards.comparison]
        ),
        json=lambda wildcards: expand(
            "fastp/json/pe/{sample}.fastp.json",
            sample=samples_per_prefixes[wildcards.comparison]
        ),
        config="multiqc/{comparison}/multiqc_config.yaml",
        fqscreen=lambda wildcards: expand(
            "fastq_screen/{sample}.{stream}.fastq_screen.{ext}",
            stream=["1", "2"],
            ext=["txt", "png"],
            sample=samples_per_prefixes[wildcards.comparison]
        ),
        additional_plots = [
            #temp("pairwise_scatterplot_mqc.png"),
            #temp("clustermap_sample_mqc.png"),
            "multiqc/{comparison}/clustermap_sample_mqc.png",
            #"multiqc/{comparison}/clustermap_genes_mqc.png",
            "multiqc/{comparison}/pca_plot_mqc.png",
            "multiqc/{comparison}/volcanoplot_mqc.png",
            "multiqc/{comparison}/distro_expr_mqc.png",
            "multiqc/{comparison}/ma_plot_mqc.png",
            "multiqc/{comparison}/distro_mu_mqc.png",
            "multiqc/{comparison}/independent_filter_mqc.png",
            "multiqc/{comparison}/inde_theta_filter_mqc.png",
            "multiqc/{comparison}/pvalue_qc_mqc.png",
            #temp("multiqc/{comparison}/clustermap_sample_mqc.png"),
            #temp("pca_axes_correlation_mqc.png")
        ]