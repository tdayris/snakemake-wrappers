.. _`deseq2_post_process`:

DESEQ2_POST_PROCESS
===================

Post process DESeq2 results and create QC report.


Example
-------

This meta-wrapper can be used by integrating the following into your workflow:

.. code-block:: python

    default_deseq2_post_process_config = {
        "condition_dict": {"DGE": {"S1": "A", "S2": "A", "S3": "B", "S4": "B"}},
        "samples_per_prefixes": {"DGE": ["S1", "S2", "S3", "S4"]},
        "genes_of_interest": ["ENSG00000141510"]
    }

    # This meta requires a list of sample per comparison level, since multiple
    # comparison may include different number of samples.
    # This dictionary should be named "samples_per_prefixes"

    # This meta requires an assignation between samples and their conditions
    # for each level of comparison. This is for colors in legends and has no
    # other impact on results.

    try:
        if config == dict():
            config = default_deseq2_post_process_config
    except NameError:
        config = default_deseq2_post_process_config


    rule all:
        input:
            reports=expand(
                "multiqc/{comparison}/MultiQC.{comparison}.html",
                comparison=config["samples_per_prefixes"].keys()
            )


    ########################
    ### Quality controls ###
    ########################

    rule multiqc:
        input:
            salmon=lambda wildcards: expand(
                "salmon/pseudo_mapping/{sample}/quant.sf",
                sample=config["samples_per_prefixes"][wildcards.comparison]
            ),
            html=lambda wildcards: expand(
                "fastp/html/pe/{sample}.fastp.html",
                sample=config["samples_per_prefixes"][wildcards.comparison]
            ),
            json=lambda wildcards: expand(
                "fastp/json/pe/{sample}.fastp.json",
                sample=config["samples_per_prefixes"][wildcards.comparison]
            ),
            config="multiqc/{comparison}/multiqc_config.yaml",
            plots = [
                #temp("pairwise_scatterplot_mqc.png"),
                #temp("clustermap_sample_mqc.png"),
                "multiqc/{comparison}/pca_plot_mqc.png",
                "multiqc/{comparison}/volcanoplot_mqc.png",
                "multiqc/{comparison}/distro_expr_mqc.png",
                "multiqc/{comparison}/distro_mu_mqc.png",
                "multiqc/{comparison}/ma_plot_mqc.png",
                "multiqc/{comparison}/independent_filter_mqc.png",
                "multiqc/{comparison}/inde_theta_filter_mqc.png",
                "multiqc/{comparison}/pvalue_qc_mqc.png"
                #temp("multiqc/{comparison}/clustermap_sample_mqc.png"),
                #temp("pca_axes_correlation_mqc.png")
            ]
        output:
            "results/{comparison}/MultiQC.{comparison}.html"
        message:
            "Aggregating quality reports from Fastp, Salmon, PCAExplorer"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: min(attempt * 1536, 10240),
            time_min=lambda wildcards, attempt: attempt * 35
        log:
            "logs/multiqc/{comparison}.log"
        params:
            lambda wildcards, input: f"--config {input.config}"
        wrapper:
            "bio/multiqc"


    rule multiqc_config:
        input:
             #pairwise_scatterplot = "image.png",
             #clustermap_sample = "image.png",
             #clustermap_genes = "figures/DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}/clustermap/ClusteredHeatmap.genes.DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}.png",
             clustermap_sample = "figures/DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}/clustermap/ClusteredHeatmap.samples.DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}.png",
             pca_plot = "figures/DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}/pca/pca_{factor}_ax_1_ax_2_with_elipse.png",
             volcanoplot = "figures/DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}/volcano/Volcano.DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}.png",
             distro_expr = "figures/DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}/log_counts/log_dst.DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}.png",
             ma_plot = "figures/DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}/maplot/maplot.DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}.png",
             distro_mu = "figures/DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}/log_counts/log_mu.DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}.png",
             independent_filter = "figures/DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}/deseq2/independent_filter.DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}.png",
             pvalue_qc = "figures/DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}/deseq2/pval.DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}.png",
             inde_theta_filter = "figures/DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}/deseq2/theta.DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}.png"
             #pca_axes_correlation = "image.png",
        output:
            multiqc_config = "multiqc/DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}/multiqc_config.yaml",
            lots = [
                #temp("pairwise_scatterplot_mqc.png"),
                #temp("multiqc/DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}/clustermap_genes_mqc.png"),
                temp("multiqc/DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}/clustermap_sample_mqc.png"),
                temp("multiqc/DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}/pca_plot_mqc.png"),
                temp("multiqc/DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}/volcanoplot_mqc.png"),
                temp("multiqc/DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}/distro_expr_mqc.png"),
                temp("multiqc/DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}/distro_mu_mqc.png"),
                temp("multiqc/DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}/ma_plot_mqc.png"),
                temp("multiqc/DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}/independent_filter_mqc.png"),
                temp("multiqc/DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}/inde_theta_filter_mqc.png"),
                temp("multiqc/DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}/pvalue_qc_mqc.png")
                #temp("multiqc/DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}/clustermap_sample_mqc.png"),
                #temp("pca_axes_correlation_mqc.png")
            ]
        message:
            "Configuring MultiQC for specialized report on {wildcards.factor} ({wildcards.test} vs {wildcards.ref})"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 512,
            time_min=lambda wildcards, attempt: attempt * 5
        log:
            "logs/multiqc/config.{factor}.{test}.{ref}.log"
        params:
            title = "Differentiel Gene Expression",
            subtitle = "Comparing {factor}: {test} VS {ref} ",
            intro_text = "This differential analysis covers {test} vs {ref}. {ref} is the reference. A fold change of 1.5 for the gene XXX means XXX is 1.5 times more expressed in {test} than in {ref}, and this difference is significative when pvalue is low (lower than 0.05).",
            report_comment = "This report has been made at Gustave Roussy.",
            show_analysis_paths = False,
            show_analysis_time = True,
            #custom_logo = '../IGR_Logo.jpeg',
            #custom_logo_url = 'https://gitlab.com/bioinfo_gustaveroussy/bigr',
            #custom_logo_title = 'BiGR, Gustave Roussy Intitute',
            report_header_info = [
                {"Contact E-mail": "bigr@gustaveroussy.fr"},
                {"Application Type": "RNA-seq"},
                {"Project Type": "Application"},
                #{"Sequencing Platform": "HiSeq 2500 High Output V4"},
                #{"Sequencing Setup": "2x125"}
            ]
        wrapper:
            "bio/BiGR/multiqc_rnaseq_report"


    ##################
    ### TSV report ###
    ##################

    rule zip_csv_report:
        input:
            "rbt/csv-report/{comparison}/html_table_deseq2_{subset}"
        output:
            report(
                "results/{comparison}/html_table_deseq2_{subset}.tar.bz2",
                caption="../report/gseapp_fc_fc.rst",
                category="9. GSEAapp Shiny",
                subcategory="{comparison}"
            )
        message:
            "Compressing {wildcards.comparison} html "
            "deseq2 {wildcards.subset} table"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 1024,
            time_min=lambda wildcards, attempt: attempt * 40,
            tmpdir="tmp"
        log:
            "logs/rbt/csv-report/compress/{comparison}_{subset}.log"
        params:
            "-cvjf"
        shell:
            "tar {params} {output} {input}"


    rule csv_report:
        input:
            "results/{comparison}/deseq2_{subset}_{comparison}.tsv"
        output:
            temp(directory("rbt/csv-report/{comparison}/html_table_deseq2_{subset}"))
        message:
            "Making {wildcards.comparison} DESeq2 results readable "
            "(DESeq2 results {wildcards.subset})"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 1024,
            time_min=lambda wildcards, attempt: attempt * 15,
            tmpdir="tmp"
        params:
            extra="--separator $'\\t' --rows-per-page 50 --sort-order ascending"
        log:
            "logs/rbt/csv-report/{comparison}_{subset}.log"
        wrapper:
            "bio/rbt/csvreport"


    rule deseq2_to_gseaapp:
        input:
            tsv = "deseq2/{comparison}/wald.{comparison}.tsv",
            gene2gene = "tximport/gene2gene.tsv"
        output:
            complete = report(
                "results/{comparison}/deseq2_complete_results_{comparison}.tsv",
                caption="../common/report/gseapp_complete.rst",
                category="6. DGE Tables",
                subcategory="{comparison}"
            ),
            fc_fc = report(
                "results/{comparison}/deseq2_sorted_on_fold_change_{comparison}.tsv",
                caption="../common/report/gseapp_fc_fc.rst",
                category="9. GSEAapp Shiny",
                subcategory="{comparison}"
            ),
            padj_fc = report(
                "results/{comparison}/deseq2_sorted_on_pval_{comparison}.tsv",
                category="9. GSEAapp Shiny",
                caption="../common/report/gseapp_padj_fc.rst",
                subcategory="{comparison}"
            )
        message:
            "Subsetting DESeq2 results for {wildcards.comparison}"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: min(attempt * 2048, 10240),
            time_min=lambda wildcards, attempt: min(attempt * 40, 200),
            tmpdir="tmp"
        params:
            ref_sample="", #TODO
            test_sample="" #TODO
        log:
            "logs/deseq2_to_gseaapp/{comparison}.log"
        wrapper:
            "bio/pandas/deseq2_to_gseaapp"



    rule deseq2_to_gseaapp_with_counts:
        input:
            tsv = "deseq2/{comparison}/wald.{comparison}.tsv",
            gene2gene = "tximport/gene2gene.tsv",
            dst = "deseq2/{comparison}/dst.{comparison}.tsv"
        output:
            complete = report(
                "results/{comparison}/deseq2_complete_results_with_counts_{comparison}.tsv",
                caption="../common/report/gseapp_complete_counts.rst",
                category="6. DGE Tables",
                subcategory="{comparison}"
            ),
            fc_fc = report(
                "results/{comparison}/deseq2_sorted_on_fold_change_with_counts_{comparison}.tsv",
                caption="../common/report/gseapp_fc_fc_counts.rst",
                category="9. GSEAapp Shiny",
                subcategory="{comparison}"
            ),
            padj_fc = report(
                "results/{comparison}/deseq2_sorted_on_pval_with_counts_{comparison}.tsv",
                category="9. GSEAapp Shiny",
                caption="../common/report/gseapp_padj_fc_counts.rst",
                subcategory="{comparison}"
            )
        message:
            "Subsetting DESeq2 results for {wildcards.comparison} with counts"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: min(attempt * 2048, 10240),
            time_min=lambda wildcards, attempt: min(attempt * 40, 200),
            tmpdir="tmp"
        log:
            "logs/deseq2_to_gseaapp/{comparison}.counts.log"
        wrapper:
            "bio/pandas/deseq2_to_gseaapp"


    ###############
    ### Seaborn ###
    ###############

    """
    This rule performs various quality control graphs and per-gene information plots
    """
    rule plot_deseq_genes:
        input:
            deseq="deseq2/{comparison}/wald.{comparison}.tsv",
            intermediar="deseq2/{comparison}/mcols.{comparison}.tsv",
            dst="deseq2/{comparison}/dst.{comparison}.tsv",
            assays="deseq2/{comparison}/assays.mu.{comparison}.tsv",
            gene2gene="tximport/gene2gene.tsv",
            metadata="deseq2/{comparison}/metadata.{comparison}.tsv",
            filter_theta="deseq2/{comparison}/filter.theta.{comparison}.tsv",
        output:
            log_counts="figures/{comparison}/log_counts/log_dst.{comparison}.png",
            log_mu="figures/{comparison}/log_counts/log_mu.{comparison}.png",
            gene_plots=report(
                expand(
                    "results/{comparison}/gene_plots/{gene}.png",
                    gene=config.get("genes_of_interest", ["ENSG00000141510"]),
                    allow_missing=True
                ),
                category="10. Gene Expression plots",
                caption="../common/reports/gene_expression_plot.rst",
                subcategory="{comparison}"
            ),
            independent_filtering="figures/{comparison}/deseq2/independent_filter.{comparison}.png",
            pval="figures/{comparison}/deseq2/pval.{comparison}.png",
            filter_theta="figures/{comparison}/deseq2/theta.{comparison}.png"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 1024,
            time_min=lambda wildcards, attempt: attempt * 15,
            tmpdir="tmp"
        params:
            condition_dict=lambda wildcards: config["condition_dict"][wildcards.comparison],
            gene_list=config.get(
                "genes_of_interest", ["ENSG00000141510"]
            ),
            gene_plots_prefix=lambda wildcards: f"results/{wildcards.comparison}/gene_plots/",
            comparison=lambda wildcards: wildcards.comparison,
            chromosomes=config.get(
                "chromosomes",
                list(range(24)) + ["MT", "X", "Y"] + list(map(str, range(24)))
            )
        log:
            "logs/deseq2/plot_genes/{comparison}.log"
        wrapper:
            "bio/seaborn/plot_deseq_genes"

    """
    This rule creates a gene-clustered heatmap
    """
    rule seaborn_clustermap_gene:
        input:
            counts = "deseq2/{comparison}/filtered_dst.{comparison}.tsv"
        output:
            png = "figures/{comparison}/clustermap/ClusteredHeatmap.genes.{comparison}.png"
        message:
            "Plotting gene-clustered heatmap for {wildcards.comparison}"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 4096,
            time_min=lambda wildcards, attempt: min(attempt * 10, 20),
            tmpdir="tmp"
        params:
            conditions=lambda wildcards: config["condition_dict"][wildcards.comparison],
            factor="{comparison}",
            row_condition_color=False
        log:
            "logs/seaborn/clustermap/{comparison}.genes.log"
        wrapper:
            "bio/seaborn/clustermap_genes"


    """
    This rule creates a sample-clustered heatmap
    """
    rule seaborn_clustermap_sample:
        input:
            counts = "deseq2/{comparison}/dst.{comparison}.tsv"
        output:
            png = "figures/{comparison}/clustermap/ClusteredHeatmap.samples.{comparison}.png"
        message:
            "Plotting sample-clustered heatmap for {wildcards.comparison}"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: min(attempt * 512, 1024),
            time_min=lambda wildcards, attempt: min(attempt * 10, 20),
            tmpdir="tmp"
        params:
            conditions=lambda wildcards: config["condition_dict"][wildcards.comparison],
            factor=lambda wildcards: (
                str(wildcards.comparison)[len("DGE_considering_factor_"):]
                if str(wildcards.comparison).startswith("DGE_considering_factor_")
                else str(wildcards.comparison)
            )
        log:
            "logs/seaborn/clustermap/{comparison}.sample.log"
        wrapper:
            "bio/seaborn/clustermap"

    rule pandas_deseq2_merge:
        input:
            wald_tsv = "deseq2/{comparison}/wald.{comparison}.tsv",
            dst_tsv = "deseq2/{comparison}/dst.{comparison}.tsv",
            gene2gene = "tximport/gene2gene.tsv"
        output:
            filtered_counts="deseq2/{comparison}/filtered_dst.{comparison}.tsv",
            filtered_deseq2="deseq2/{comparison}/filtered_deseq2.{comparison}.tsv",
            merged_table="deseq2/{comparison}/merged_deseq2_counts.{comparison}.tsv"
        message:
            "Merging counts and DESeq2 results for {wildcards.comparison}"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 1024 * 6,
            time_min=lambda wildcards, attempt: min(attempt * 10, 20),
            tmpdir="tmp"
        params:
            alpha_threshold=0.05,
            fold_change=0.001
        log:
            "logs/pandas/merged_deseq2_counts/{comparison}.log"
        wrapper:
            "bio/pandas/deseq2_merge"

    #######################
    ### EnhancedVolcano ###
    #######################

    """
    This rules computes and plots a Volcano-plot
    """
    rule enhancedvolcano_volcanoplot:
        input:
            deseq2_tsv="deseq2/{comparison}/wald.{comparison}.tsv"
        output:
            png="figures/{comparison}/volcano/Volcano.{comparison}.png"
        message: "Plotting Volcanoplot for {wildcards.comparison}"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 2048,
            time_min=lambda wildcards, attempt: attempt * 15,
            tmpdir="tmp"
        params:
            alpha_threshold=config["thresholds"].get("alpha", 0.05),
            fc_threshold=config["thresholds"].get("fc", 0.6)
        log:
            "logs/enhanced_volcano/{comparison}.log"
        wrapper:
            "bio/enhancedVolcano/volcano-deseq2"


    ############################
    ### Consensus Clustering ###
    ############################

    rule consensus_cluster_plus:
        input:
            expr_mat = "deseq2/{comparison}/dst.{comparison}.tsv"
        output:
            res_dir=directory("consensusclusterplus/{comparison}")
        message:
            "Clustering samples in {wildcards.comparison}"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 1024,
            time_min=lambda wildcards, attempt: attempt * 15,
            tmpdir="tmp"
        log:
            "logs/consensus_cluster_plus/{comparison}.log"
        params:
            extra=lambda wildcards: f"clusterAlg='km', title='{wildcards.comparison}', distance='euclidean'"
        wrapper:
            "bio/consensusclusterplus"


    #######################################
    ### General PCA over all the cohort ###
    #######################################

    """
    This rule plots a PCA on a given TSV count file
    """
    rule general_pca:
        input:
            counts="salmon/TMP.genes.nochr.tsv"
        output:
            png=expand(
                "figures/pca/general.pca_{axes}.png",
                axes=["PC1_PC2", "PC2_PC1"]
            )
        message:
            "Plotting general PCA"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 4096,
            time_min=lambda wildcards, attempt: attempt * 5,
            tmpdir="tmp"
        log:
            "logs/seaborn/pca/general.log"
        params:
            axes=[1, 2],
            # conditions=lambda wildcards: dict(zip(
            #     config["design"].index.tolist(),
            #     config["design"][wildcards.factor].tolist()
            # )),
            conditions = dict(zip(
                config["design"].index.tolist(),
                config["design"].drop(["Sample_id", "Upstream_file", "Downstream_file"], axis=1).iloc[:, :-1].apply("_".join, axis=1).to_list()
            )),
            prefix=lambda wildcards: f"figures/pca/general.pca"
        wrapper:
            "bio/seaborn/pca"


    """
    This rule removes unnecessary columns for further seaborn PCA
    """
    rule filter_merged_counts:
        input:
            table="salmon/TPM.genes.tsv"
        output:
            table=temp("salmon/TMP.genes.nochr.tsv")
        message:
            "Removing gene annotations"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 4096,
            time_min=lambda wildcards, attempt: attempt * 5,
            tmpdir="tmp"
        group:
            "salmon_merge"
        log:
            "logs/seaborn/pca/filter.log"
        params:
            drop_column=["Chromosomes", "Start", "End", "Hugo_id"]
        wrapper:
            "bio/pandas/filter_table"


    """
    This rule merges salmon counts into a single TSV table
    """
    rule pandas_merge_salmon_tr:
        input:
            quant = expand(
                "salmon/pseudo_mapping/{sample}/quant.sf",
                sample=config["design"].index.tolist()
            ),
            tx2gene = "tximport/transcripts2genes.tsv"
        output:
            tsv = "salmon/TPM.{content}.tsv"
        message:
            "Merging salmon individual counts on {wildcards.content} "
            "for further general PCA"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 3072,
            time_min=lambda wildcards, attempt: attempt * 5,
            tmpdir="tmp"
        group:
            "salmon_merge"
        params:
            header = True,
            position = True,
            gencode = True,
            drop_na = True,
            dro_null = True,
            genes = lambda wildcards: wildcards.content == "genes"
        log:
            "logs/pandas_merge_salmon/{content}.log"
        wrapper:
            "bio/pandas/salmon"


    ####################
    ### PCA Explorer ###
    ####################

    """
    This rule simply plots the PCA
    """
    rule pcaexplorer_pca:
        input:
            dst = "deseq2/DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}/wald.DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}.RDS"
        output:
            png = "figures/DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}/pca/pca_{factor}_ax_{a}_ax_{b}_{elipse}.png"
        message:
            "Plotting PCA for ({wildcards.factor}:"
            "{wildcards.a}/{wildcards.b}:{wildcards.elipse})"
        threads:
            1
        resources:
            mem_mb=lambda wildcards, attempt: min(attempt * 1024, 10240),
            time_min=lambda wildcards, attempt: min(attempt * 20, 200),
            tmpdir="tmp"
        params:
            extra=lambda wildcards: (
                f"intgroup = c('{wildcards.factor}'), ntop = 100, pcX = {wildcards.a}, pcY = {wildcards.b}, ellipse = {'TRUE' if wildcards.elipse == 'with_elipse' else 'FALSE'}"),
            w = 1024,
            h = 768
        log:
            "logs/pcaexplorer/PCA/DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}/pca_ingroup_{factor}_ax_{a}_{b}_{elipse}.log"
        wrapper:
            "bio/pcaExplorer/PCA"


    """
    This rule plots the distribution of the expression of Salmon counts
    """
    rule pca_explorer_distro_expr:
        input:
            dst = "deseq2/{comparison}/wald.{comparison}.RDS"
        output:
            png = "figures/{comparison}/distro_expr/distro_expr.{comparison}.png"
        message:
            "Plotting expression distributions for {wildcards.comparison}"
        threads: 1
        resources:
            mem_mb = (
                lambda wildcards, attempt: min(attempt * 1024, 10240)
            ),
            time_min = (
                lambda wildcards, attempt: min(attempt * 20, 200)
            )
        log:
            "logs/pcaexplorer/distro_expr/{comparison}.log"
        wrapper:
            "bio/pcaExplorer/distro_expr"


    ##############
    ### DESeq2 ###
    ##############
    """
    This rule creates a MA-Plot
    """
    rule deseq2_maplot:
        input:
            res = "deseq2/{comparison}/wald.{comparison}.tsv"
        output:
            png = "figures/{comparison}/maplot/maplot.{comparison}.png"
        message:
            "Building MA-plot for {wildcards.comparison}"
        threads: 1
        resources:
            mem_mb = (
                lambda wildcards, attempt: min(attempt * 1024, 10240)
            ),
            time_min = (
                lambda wildcards, attempt: min(attempt * 20, 200)
            )
        log:
            "logs/deseq2/maplot/maplot.{comparison}.log"
        wrapper:
            "bio/deseq2/plotMA"

Note that input, output and log file paths can be chosen freely, as long as the dependencies between the rules remain as listed here.
For additional parameters in each individual wrapper, please refer to their corresponding documentation (see links below).

When running with

.. code-block:: bash

    snakemake --use-conda

the software dependencies will be automatically deployed into an isolated environment before execution.



Used wrappers
---------------------

The following individual wrappers are used in this meta-wrapper:


* :ref:`bio/multiqc`

* :ref:`bio/BiGR/multiqc_rnaseq_report`

* :ref:`bio/rbt/csvreport`

* :ref:`bio/pandas/deseq2_to_gseaapp`

* :ref:`bio/seaborn/plot_deseq_genes`

* :ref:`bio/seaborn/clustermap_genes`

* :ref:`bio/seaborn/clustermap`

* :ref:`bio/pandas/deseq2_merge`

* :ref:`bio/consensusclusterplus`

* :ref:`bio/enhancedVolcano/volcano-deseq2`

* :ref:`bio/seaborn/pca`

* :ref:`bio/pandas/filter_table`

* :ref:`bio/pandas/salmon`

* :ref:`bio/pcaExplorer/PCA`

* :ref:`bio/pcaExplorer/distro_expr`

* :ref:`bio/deseq2/plotMA`


Please refer to each wrapper in above list for additional configuration parameters and information about the executed code.







Authors
-------


* Thibault Dayris

