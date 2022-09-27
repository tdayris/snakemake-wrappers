##############
### Tables ###
##############

# This rule makes the DEseq2 results human-readable
# (and compatible with gseaapp from BiGR)
"""
009.deseq2_readable
from:
-> 004.tx_to_gene
-> 008.deseq2
by:
->
"""
rule 009_deseq2_readable:
    input:
        tsv="008.deseq2/{comparison}/wald.{comparison}.tsv",
        gene2gene="004.salmon/gene2gene_with_chr.tsv",
        dst="008.deseq2/{comparison}/dst.{comparison}.tsv",
    output:
        complete=report(
            "data_output/DEseq2/{comparison}/Complete_{comparison}.tsv",
            caption=str(workflow_source_dir / "reports" / "009.gseapp_complete.rst"),
            category="DEseq2",
            subcategory="{comparison}",
        ),
        fc_fc=report(
            "data_output/DEseq2/{comparison}/SortedOnLogFC_{comparison}.tsv",
            caption=str(workflow_source_dir / "reports" / "009.gseapp_fc_fc.rst"),
            category="DEseq2",
            subcategory="{comparison}",
        ),
        padj_fc=report(
            "data_output/DEseq2/{comparison}/SortedOnPadj_{comparison}.tsv",
            category="DESeq2",
            caption=str(workflow_source_dir / "reports" / "009.gseapp_padj_fc.rst"),
            subcategory="{comparison}",
        ),
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_20min_per_attempt,
        tmpdir="tmp",
    retries: 1
    log:
        "logs/009.deseq2/readable/{comparison}.log",
    wrapper:
        "bio/pandas/deseq2_to_gseaapp"


###################
### CSV reports ###
###################

# This rule creates HTML reports from DESeq2 results
"""
009.rbt_csv_report
from:
-> 009.deseq2_readable
by:
-> 009.zip_csv_report
"""
rule 009_rbt_csv_report:
    input:
        "data_output/DEseq2/{comparison}/{content}_{comparison}.tsv",
    output:
        temp(directory("009.rbt/csvreport/{comparison}/{content}_html_table/")),
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_15min_per_attempt,
        tmpdir="tmp",
    retries: 1
    group:
        "csv_report"
    log:
        "logs/009.rbt/csv-report/{comparison}.{content}.log",
    params:
        config["rbt"].get(
            "csv_extra",
        '--separator "\\t" --sort-column stat_change '
            "--sort-order ascending --rows-per-page 50",
        ),
    wrapper:
        "bio/rbt/csvreport"

"""
009.zip_csv_report
from:
-> 009.rbt_csv_report
by:
-> End job
"""
rule 009_zip_csv_report:
    input:
        "009.rbt/csvreport/{comparison}/{content}_html_table/",
    output:
        protected("data_output/DEseq2/{comparison}/{content}_html_table.tar.bz2"),
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir="tmp",
    retries: 1
    group:
        "csv_report"
    log:
        "logs/009.rbt/csv-report/compress/{comparison}.{content}.log",
    params:
        config["rbt"].get("zip_extra", "--create --bzip2"),
    conda:
        str(workflow_source_dir / "envs" / "bash.yaml")
    shell:
        "tar {params} --file={output} {input} > {log} 2>&1 "


################
### MultiQC  ###
################
"""
009.multiqc_config
from:
-> 009.plot_deseq_genes
by:
-> Snakefile.deseq2_results
"""
rule 009_multiqc_config:
    input:
        clustermap_sample="009.figures/DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}/clustermap/ClusteredHeatmap.samples.DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}.png",
        pca_plot="009.figures/DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}/pca/pca_{factor}_ax_1_ax_2_with_elipse.png",
        volcanoplot="009.figures/DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}/volcano/Volcano.DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}.png",
        distro_expr="009.figures/DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}/log_counts/log_dst.DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}.png",
        ma_plot="009.figures/DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}/maplot/maplot.DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}.png",
        distro_mu="009.figures/DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}/log_counts/log_mu.DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}.png",
        independent_filter="009.figures/DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}/deseq2/independent_filter.DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}.png",
        pvalue_qc="009.figures/DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}/deseq2/pval.DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}.png",
        inde_theta_filter="009.figures/DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}/deseq2/theta.DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}.png",
    output:
        multiqc_config=temp(
            "009.multiqc/DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}/multiqc_config.yaml"
        ),
        lots=[
            temp(
                "009.multiqc/DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}/clustermap_sample_mqc.png"
            ),
            temp(
                "009.multiqc/DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}/pca_plot_mqc.png"
            ),
            temp(
                "009.multiqc/DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}/volcanoplot_mqc.png"
            ),
            temp(
                "009.multiqc/DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}/distro_expr_mqc.png"
            ),
            temp(
                "009.multiqc/DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}/distro_mu_mqc.png"
            ),
            temp(
                "009.multiqc/DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}/ma_plot_mqc.png"
            ),
            temp(
                "009.multiqc/DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}/independent_filter_mqc.png"
            ),
            temp(
                "009.multiqc/DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}/inde_theta_filter_mqc.png"
            ),
            temp(
                "009.multiqc/DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}/pvalue_qc_mqc.png"
            ),
        ],
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_15min_per_attempt,
        tmpdir="tmp",
    retries: 1
    log:
        "logs/009.multiqc/config.{factor}.{test}.{ref}.log",
    params:
        title=config["multiqc"].get("title", "Differentiel Gene Expression"),
        subtitle="Comparing {factor}: {test} VS {ref} ",
        intro_text="This differential analysis covers {test} vs {ref}. {ref} is the reference. A fold change of 1.5 for the gene XXX means XXX is 1.5 times more expressed in {test} than in {ref}, and this difference is significative when pvalue is low (lower than 0.05).",
        report_comment=config["multiqc"].get(
            "report_comment", "This report has been made at Gustave Roussy."
        ),
        show_analysis_paths=False,
        show_analysis_time=True,
        report_header_info=[
            {
                "Contact E-mail": config["multiqc"].get(
        "email", "bigr@gustaveroussy.fr"
                )
            },
            {"Application Type": config["multiqc"].get("application_type", "RNA-seq")},
            {"Project Type": config["multiqc"].get("project_type", "Application")},
        ],
    wrapper:
        "bio/BiGR/multiqc_rnaseq_report"


###############
### Seaborn ###
###############

# This rule performs various quality control graphs and per-gene information plots
"""

"""
rule 009_plot_deseq_genes:
    input:
        deseq="008.deseq2/{comparison}/wald.{comparison}.tsv",
        intermediar="008.deseq2/{comparison}/mcols.{comparison}.tsv",
        dst="008.deseq2/{comparison}/dst.{comparison}.tsv",
        assays="008.deseq2/{comparison}/assays.mu.{comparison}.tsv",
        gene2gene="004.salmon/gene2gene_with_chr.tsv",
        metadata="008.deseq2/{comparison}/metadata.{comparison}.tsv",
        filter_theta="008.deseq2/{comparison}/filter.theta.{comparison}.tsv",
    output:
        log_counts=temp("009.figures/{comparison}/log_counts/log_dst.{comparison}.png"),
        log_mu=temp("009.figures/{comparison}/log_counts/log_mu.{comparison}.png"),
        gene_plots=directory("data_output/DEseq2/{comparison}/gene_plots/"),
        independent_filtering=temp(
            "009.figures/{comparison}/deseq2/independent_filter.{comparison}.png"
        ),
        pval=temp("009.figures/{comparison}/deseq2/pval.{comparison}.png"),
        filter_theta=temp("009.figures/{comparison}/deseq2/theta.{comparison}.png"),
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_15min_per_attempt,
        tmpdir="tmp",
    retries: 1
    params:
        condition_dict=lambda wildcards: condition_dict[wildcards.comparison],
        gene_list=config.get("genes_of_interest", ["ENSG00000141510"]),
        gene_plots_prefix=lambda wildcards: f"data_output/DEseq2/{wildcards.comparison}/gene_plots/",
        comparison=lambda wildcards: wildcards.comparison,
        chromosomes=config["reference"].get(
            "chromosomes",
            list(range(24)) + ["MT", "X", "Y"] + list(map(str, range(24))),
        ),
    log:
        "logs/009.deseq2/plot_genes/{comparison}.log",
    wrapper:
        "bio/seaborn/plot_deseq_genes"


"""
This rule creates a sample-clustered heatmap
"""


rule 009_seaborn_clustermap_sample:
    input:
        counts="008.deseq2/{comparison}/dst.{comparison}.tsv",
    output:
        png=temp(
            "009.figures/{comparison}/clustermap/ClusteredHeatmap.samples.{comparison}.png"
        ),
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_15min_per_attempt,
        tmpdir="tmp",
    retries: 1
    params:
        conditions=lambda wildcards: condition_dict[wildcards.comparison],
        factor=lambda wildcards: (
            str(wildcards.comparison)[len("DGE_considering_factor_") :]
            if str(wildcards.comparison).startswith("DGE_considering_factor_")
            else str(wildcards.comparison)
        ),
    log:
        "logs/009.seaborn/clustermap/{comparison}.sample.log",
    wrapper:
        "bio/seaborn/clustermap"


#######################
### EnhancedVolcano ###
#######################

"""
This rules computes and plots a Volcano-plot
"""


rule 009_enhancedvolcano_volcanoplot:
    input:
        deseq2_tsv="008.deseq2/{comparison}/wald.{comparison}.tsv",
    output:
        png=temp("009.figures/{comparison}/volcano/Volcano.{comparison}.png"),
    threads: 1
    resources:
        mem_mb=get_2gb_per_attempt,
        time_min=get_15min_per_attempt,
        tmpdir="tmp",
    retries: 1
    params:
        alpha_threshold=config["deseq2"]["thresholds"].get("alpha", 0.05),
        fc_threshold=config["deseq2"]["thresholds"].get("fc", 0.6),
    log:
        "logs/009.enhanced_volcano/{comparison}.log",
    wrapper:
        "bio/enhancedVolcano/volcano-deseq2"


##############
### DESeq2 ###
##############
"""
This rule creates a MA-Plot
"""


rule 009_deseq2_maplot:
    input:
        res="008.deseq2/{comparison}/wald.{comparison}.tsv",
    output:
        png=temp("009.figures/{comparison}/maplot/maplot.{comparison}.png"),
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_15min_per_attempt,
        tmpdir="tmp",
    retries: 1
    log:
        "logs/009.deseq2/maplot/maplot.{comparison}.log",
    wrapper:
        "bio/deseq2/plotMA"


####################
### PCA Explorer ###
####################

"""
This rule simply plots the PCA
"""


rule 009_pcaexplorer_pca:
    input:
        dst="008.deseq2/DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}/wald.DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}.RDS",
    output:
        png=temp(
            "009.figures/DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}/pca/pca_{factor}_ax_{a}_ax_{b}_{elipse}.png"
        ),
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_15min_per_attempt,
        tmpdir="tmp",
    retries: 1
    params:
        extra=lambda wildcards: (
            f"intgroup = c('{wildcards.factor}'), ntop = 100, pcX = {wildcards.a}, pcY = {wildcards.b}, ellipse = {'TRUE' if wildcards.elipse == 'with_elipse' else 'FALSE'}"
        ),
        w=1024,
        h=768,
    log:
        "logs/009.pcaexplorer/PCA/DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}/pca_ingroup_{factor}_ax_{a}_{b}_{elipse}.log",
    wrapper:
        "bio/pcaExplorer/PCA"


"""
This rule plots the distribution of the expression of Salmon counts
"""


rule 009_pca_explorer_distro_expr:
    input:
        dst="008.deseq2/{comparison}/wald.{comparison}.RDS",
    output:
        png=temp("009.figures/{comparison}/distro_expr/distro_expr.{comparison}.png"),
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_15min_per_attempt,
        tmpdir="tmp",
    retries: 1
    log:
        "logs/009.pcaexplorer/distro_expr/{comparison}.log",
    wrapper:
        "bio/pcaExplorer/distro_expr"
