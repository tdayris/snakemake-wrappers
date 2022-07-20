##############
### Tables ###
##############
"""
This rule makes the DEseq2 results human-readable
(and compatible with gseaapp from BiGR)
"""


rule deseq2_readable:
    input:
        tsv="deseq2/{comparison}/wald.{comparison}.tsv",
        gene2gene="tximport/gene2gene.tsv",
        dst="deseq2/{comparison}/dst.{comparison}.tsv",
    output:
        complete=report(
            "data_output/DEseq2/{comparison}/Complete_{comparison}.tsv",
            caption=str(worflow_source_dir / "reports" / "009.gseapp_complete.rst"),
            category="DEseq2",
            subcategory="{comparison}",
        ),
        fc_fc=report(
            "data_output/DEseq2/{comparison}/SortedOnLogFC_{comparison}.tsv",
            caption=str(worflow_source_dir / "reports" / "009.gseapp_fc_fc.rst"),
            category="DEseq2",
            subcategory="{comparison}",
        ),
        padj_fc=report(
            "data_output/DEseq2/{comparison}/SortedOnPadj_{comparison}.tsv",
            category="DESeq2",
            caption=str(worflow_source_dir / "reports" / "009.gseapp_padj_fc.rst"),
            subcategory="{comparison}",
        ),
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_20min_per_attempt,
        tmpdir="tmp",
    log:
        "logs/deseq2/readable/{comparison}.log",
    wrapper:
        "bio/pandas/deseq2_to_gseaapp"


###################
### CSV reports ###
###################

"""
This rule creates HTML reports from DESeq2 results
"""


rule rbt_csv_report:
    input:
        "deseq2/{comparison}/{comparison}.{content}.tsv",
    output:
        directory("data_output/DEseq2/{comparison}/{comparison}.{content}/"),
    message:
        "Making DESeq2 results readable and searchable"
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_15min_per_attempt,
        tmpdir="tmp",
    group:
        "csv_report"
    log:
        "logs/rbt/csv-report/{comparison}.{content}.log",
    params:
        '--separator "\\t" '
        "--sort-column stat_change "
        "--sort-order ascending "
        "--rows-per-page 50",
    wrapper:
        "bio/rbt/csvreport"


rule zip_csv_report:
    input:
        "rbt/csv-report/{comparison}/html_table_deseq2_{subset}",
    output:
        "data_output/DESeq2/{comparison}/html_table_deseq2_{subset}.tar.bz2",
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir="tmp",
    group:
        "csv_report"
    log:
        "logs/rbt/csv-report/compress/{comparison}_{subset}.log",
    params:
        "-cvjf",
    shell:
        "tar {params} {output} {input}"


################
### MultiQC  ###
################


rule multiqc_config:
    input:
        clustermap_sample="figures/DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}/clustermap/ClusteredHeatmap.samples.DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}.png",
        pca_plot="figures/DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}/pca/pca_{factor}_ax_1_ax_2_with_elipse.png",
        volcanoplot="figures/DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}/volcano/Volcano.DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}.png",
        distro_expr="figures/DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}/log_counts/log_dst.DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}.png",
        ma_plot="figures/DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}/maplot/maplot.DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}.png",
        distro_mu="figures/DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}/log_counts/log_mu.DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}.png",
        independent_filter="figures/DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}/deseq2/independent_filter.DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}.png",
        pvalue_qc="figures/DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}/deseq2/pval.DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}.png",
        inde_theta_filter="figures/DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}/deseq2/theta.DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}.png",
    output:
        multiqc_config=temp(
            "multiqc/DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}/multiqc_config.yaml"
        ),
        lots=[
            temp(
                "multiqc/DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}/clustermap_sample_mqc.png"
            ),
            temp(
                "multiqc/DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}/pca_plot_mqc.png"
            ),
            temp(
                "multiqc/DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}/volcanoplot_mqc.png"
            ),
            temp(
                "multiqc/DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}/distro_expr_mqc.png"
            ),
            temp(
                "multiqc/DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}/distro_mu_mqc.png"
            ),
            temp(
                "multiqc/DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}/ma_plot_mqc.png"
            ),
            temp(
                "multiqc/DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}/independent_filter_mqc.png"
            ),
            temp(
                "multiqc/DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}/inde_theta_filter_mqc.png"
            ),
            temp(
                "multiqc/DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}/pvalue_qc_mqc.png"
            ),
        ],
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_15min_per_attempt,
        tmpdir="tmp",
    log:
        "logs/multiqc/config.{factor}.{test}.{ref}.log",
    params:
        title="Differentiel Gene Expression",
        subtitle="Comparing {factor}: {test} VS {ref} ",
        intro_text="This differential analysis covers {test} vs {ref}. {ref} is the reference. A fold change of 1.5 for the gene XXX means XXX is 1.5 times more expressed in {test} than in {ref}, and this difference is significative when pvalue is low (lower than 0.05).",
        report_comment="This report has been made at Gustave Roussy.",
        show_analysis_paths=False,
        show_analysis_time=True,
        report_header_info=[
            {"Contact E-mail": "bigr@gustaveroussy.fr"},
            {"Application Type": "RNA-seq"},
            {"Project Type": "Application"},
        ],
    wrapper:
        "bio/BiGR/multiqc_rnaseq_report"


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
            allow_missing=True,
        ),
        category="Gene Expression plots",
        caption=str(
        worflow_source_dir / "reports" / " 009.gene_expression_plot.rst"
            ),
            subcategory="{comparison}",
        ),
        independent_filtering="figures/{comparison}/deseq2/independent_filter.{comparison}.png",
        pval="figures/{comparison}/deseq2/pval.{comparison}.png",
        filter_theta="figures/{comparison}/deseq2/theta.{comparison}.png",
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_15min_per_attempt,
        tmpdir="tmp",
    params:
        condition_dict=lambda wildcards: condition_dict[wildcards.comparison],
        gene_list=config.get("genes_of_interest", ["ENSG00000141510"]),
        gene_plots_prefix=lambda wildcards: f"results/{wildcards.comparison}/gene_plots/",
        comparison=lambda wildcards: wildcards.comparison,
        chromosomes=config["reference"].get(
            "chromosomes",
            list(range(24)) + ["MT", "X", "Y"] + list(map(str, range(24))),
        ),
    log:
        "logs/deseq2/plot_genes/{comparison}.log",
    wrapper:
        "bio/seaborn/plot_deseq_genes"


"""
This rule creates a sample-clustered heatmap
"""


rule seaborn_clustermap_sample:
    input:
        counts="deseq2/{comparison}/dst.{comparison}.tsv",
    output:
        png="figures/{comparison}/clustermap/ClusteredHeatmap.samples.{comparison}.png",
    message:
        "Plotting sample-clustered heatmap for {wildcards.comparison}"
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_15min_per_attempt,
        tmpdir="tmp",
    params:
        conditions=lambda wildcards: condition_dict[wildcards.comparison],
        factor=lambda wildcards: (
            str(wildcards.comparison)[len("DGE_considering_factor_") :]
            if str(wildcards.comparison).startswith("DGE_considering_factor_")
            else str(wildcards.comparison)
        ),
    log:
        "logs/seaborn/clustermap/{comparison}.sample.log",
    wrapper:
        "bio/seaborn/clustermap"


#######################
### EnhancedVolcano ###
#######################

"""
This rules computes and plots a Volcano-plot
"""


rule enhancedvolcano_volcanoplot:
    input:
        deseq2_tsv="deseq2/{comparison}/wald.{comparison}.tsv",
    output:
        png="figures/{comparison}/volcano/Volcano.{comparison}.png",
    message:
        "Plotting Volcanoplot for {wildcards.comparison}"
    threads: 1
    resources:
        mem_mb=get_2gb_per_attempt,
        time_min=get_15min_per_attempt,
        tmpdir="tmp",
    params:
        alpha_threshold=config["deseq2"]["thresholds"].get("alpha", 0.05),
        fc_threshold=config["deseq2"]["thresholds"].get("fc", 0.6),
    log:
        "logs/enhanced_volcano/{comparison}.log",
    wrapper:
        "bio/enhancedVolcano/volcano-deseq2"


##############
### DESeq2 ###
##############
"""
This rule creates a MA-Plot
"""


rule deseq2_maplot:
    input:
        res="deseq2/{comparison}/wald.{comparison}.tsv",
    output:
        png="figures/{comparison}/maplot/maplot.{comparison}.png",
    message:
        "Building MA-plot for {wildcards.comparison}"
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_15min_per_attempt,
        tmpdir="tmp",
    log:
        "logs/deseq2/maplot/maplot.{comparison}.log",
    wrapper:
        "bio/deseq2/plotMA"


####################
### PCA Explorer ###
####################

"""
This rule simply plots the PCA
"""


rule pcaexplorer_pca:
    input:
        dst="deseq2/DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}/wald.DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}.RDS",
    output:
        png="figures/DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}/pca/pca_{factor}_ax_{a}_ax_{b}_{elipse}.png",
    message:
        "Plotting PCA for ({wildcards.factor}:"
        "{wildcards.a}/{wildcards.b}:{wildcards.elipse})"
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_15min_per_attempt,
        tmpdir="tmp",
    params:
        extra=lambda wildcards: (
            f"intgroup = c('{wildcards.factor}'), ntop = 100, pcX = {wildcards.a}, pcY = {wildcards.b}, ellipse = {'TRUE' if wildcards.elipse == 'with_elipse' else 'FALSE'}"
        ),
        w=1024,
        h=768,
    log:
        "logs/pcaexplorer/PCA/DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}/pca_ingroup_{factor}_ax_{a}_{b}_{elipse}.log",
    wrapper:
        "bio/pcaExplorer/PCA"


"""
This rule plots the distribution of the expression of Salmon counts
"""


rule pca_explorer_distro_expr:
    input:
        dst="deseq2/{comparison}/wald.{comparison}.RDS",
    output:
        png="figures/{comparison}/distro_expr/distro_expr.{comparison}.png",
    message:
        "Plotting expression distributions for {wildcards.comparison}"
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_15min_per_attempt,
        tmpdir="tmp",
    log:
        "logs/pcaexplorer/distro_expr/{comparison}.log",
    wrapper:
        "bio/pcaExplorer/distro_expr"
