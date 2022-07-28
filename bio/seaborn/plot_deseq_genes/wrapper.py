#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Draw several QC plots for DESeq2
"""

import logging
import matplotlib
import matplotlib.pyplot
import numpy
import pandas
import seaborn

from typing import Any, Dict, List, Optional

logging.basicConfig(
    #filename=snakemake.log[0],
    filemode="w",
    level=logging.DEBUG
)

logging.getLogger('matplotlib.font_manager').disabled = True

def read_deseq(path: str) -> pandas.DataFrame:
    logging.info(f"Loading {path}")
    tmp = pandas.read_csv(
        path,
        sep="\t",
        header=0,
        index_col=0
    )
    tmp.index.name = "Ensembl_ID"
    return tmp

def read_gene2gene(path: str) -> pandas.DataFrame:
    logging.info(f"Loading {path}")
    tmp = pandas.read_csv(
        path,
        sep="\t",
        header=0,
        index_col=3
    )
    tmp.index.name = "Ensembl_ID"
    return tmp


def load_all_data(deseq_path: str,
                  intermediar_path: str,
                  assays_path: str,
                  dst_path: str,
                  filter_theta_path: str,
                  metadata_path: str,
                  gene2gene_path: str,
                  padj_threshold: float = 0.05) -> List[pandas.DataFrame]:
    deseq = pandas.merge(
        left=read_deseq(deseq_path),
        right=read_deseq(intermediar_path),
        left_index=True,
        right_index=True,
        how="outer",
        validate="1:1",
        suffixes=["_results", "_mcols"]
    )
    deseq = pandas.merge(
        left=deseq,
        right=read_gene2gene(gene2gene_path),
        left_index=True,
        right_index=True,
        how="outer",
        validate="1:1",
    )
    deseq["Status"] = [
        "Differentially Expressed"
        if padj <= padj_threshold else "Not Significative"
        for padj in deseq["padj"]
    ]
    deseq["log2FC_eq_lfcSE"] = [
        (fc > 0) is (shr > 0)
        for fc, shr in zip(deseq["log2FoldChange"], deseq["lfcSE"])
    ]
    print(deseq.columns)

    counts = read_deseq(dst_path)
    mu = read_deseq(assays_path)
    metadata = read_deseq(metadata_path)
    theta = read_deseq(filter_theta_path)
    return deseq, counts, mu, theta, metadata


def plot_counts(dst: pandas.DataFrame,
                colors: Dict[str, str],
                png_out: str,
                is_mu: bool = False) -> None:
    """
    Boxplot of log10 counts of each sample,
    colored by conditions
    """
    count = "Counts" if is_mu is False else "µ"
    logging.info(f"Plotting log10 {count}")
    box = dst.loc[~(dst==0).all(axis=1)]
    box = box.melt(
        ignore_index=False,
        # Counts are NOT log10, but they will be in final
        # graph !
        value_name=f"log10(Normalised {count})",
        var_name="Sample ID"
    )

    # g, ax = matplotlib.pyplot.subplots(1, 1)
    g = seaborn.catplot(
        # Counts are NOT log10, but they will be in final
        # graph !
        x=f"log10(Normalised {count})",
        y="Sample ID",
        data=box,
        palette=colors,
        kind="box",
        legend=True,
        legend_out=True
    )
    g.set(xscale="log")

    matplotlib.pyplot.savefig(
        png_out,
        bbox_inches="tight"
    )


def plot_single_gene(dst: pandas.DataFrame,
                     mu: pandas.DataFrame,
                     deseq: pandas.DataFrame,
                     condition_dict: Dict[str, str],
                     genes: List[str],
                     png_prefix: str,
                     comparison: str) -> None:
    """
    Plot gene expression, its geometric mean and dispersion value
    """
    dstbox = dst.melt(
        ignore_index=False,
        value_name="Normalised Counts",
        var_name="Sample ID"
    )
    dstbox["Condition"] = [condition_dict[s] for s in dstbox['Sample ID']]

    mubox = mu.melt(
        ignore_index=False,
        value_name="µ",
        var_name="Sample ID"
    )
    mubox["Condition"] = [condition_dict[s] for s in mubox['Sample ID']]

    for gene in genes:
        logging.info(f"Plotting information on {gene}")
        gene_count = dstbox[dstbox.index == gene]
        gene_mu = mubox[mubox.index == gene]
        gene_dge = deseq[deseq.index == gene]
        png_out = png_prefix + f"{gene}.png"
        status = gene_dge["Status"][0]

        padj = round(gene_dge["padj"][0], 4)

        try:
            name = gene_dge["Gene_Name"][0]
            name = f"{name} ({gene})"
        except KeyError:
            name = gene

        fig = matplotlib.pyplot.figure(figsize=(20,10))

        axes = []
        axes.append(fig.add_subplot(1, 4, 1))
        axes.append(fig.add_subplot(1, 4, 2, sharey=axes[0]))
        axes.append(fig.add_subplot(1, 4, 3, sharey=axes[0]))
        axes.append(fig.add_subplot(1, 4, 4))

        # fig, axes = matplotlib.pyplot.subplots(
        #     1, 4, figsize=(20,10)#, sharey=True
        # )
        fig.suptitle(
            f"Expression of {name} among samples\n"
            "Differential gene expression on "
            f"{comparison} {status} with P-Adj of {padj}"
        )
        axes[0].set_title("Gene VST accross conditions")
        axes[1].set_title("Gene µ accross conditions")
        axes[2].set_title("Geometric Mean and dispersion")
        axes[3].set_title("Expression change with standard Error")
        seaborn.boxplot(
            ax=axes[0],
            y="Normalised Counts",
            x="Condition",
            data=gene_count
        )
        seaborn.boxplot(
            ax=axes[1],
            y="µ",
            x="Condition",
            data=gene_mu
        )
        axes[2].bar(
            x=[1],
            height=gene_dge["baseMean_results"][0],
            yerr=gene_dge["baseMean_results"][0] * gene_dge["dispersion"][0],
            tick_label=gene
        )
        axes[3].bar(
            x=[1],
            height=gene_dge["log2FoldChange"][0],
            tick_label="Log2(FoldChange)",
            yerr=gene_dge["lfcSE"][0]
        )

        # Rotate condition texts
        for ax in axes:
            matplotlib.pyplot.sca(ax)
            matplotlib.pyplot.xticks(rotation=90)

        matplotlib.pyplot.savefig(
            png_out,
            bbox_inches="tight"
        )


def independent_filtering(deseq: pandas.DataFrame,
                          comparison: str,
                          png_out: str) -> None:
    """Plot number of observations filtered by independent filtering method"""
    tmp = deseq[
        ["padj", "Status", "filterThreshold", "baseMean_results", "maxCooks"]
    ]
    tmp["padj"] = -numpy.log10(tmp["padj"])
    tmp["log_baseMean_results"] = numpy.log10(tmp["baseMean_results"])
    tmp["filterThreshold"] = [
        "Pass" if filter is True else "Do not pass"
        for filter in tmp["filterThreshold"]
    ]
    tmp.columns = [
        "-Log10(P-Adjusted)",
        "Status",
        "Independent Filtering",
        "Mean Expression",
        "Maximum Cook's distance",
        "Log10(Mean Expression)",
    ]

    # Build plot
    fig, axes = matplotlib.pyplot.subplots(
        1, 3, figsize=(20,10)#, sharey=True
    )
    fig.suptitle(f"Independent filtering, {comparison}")
    axes[0].set_title("# observations filtered out")
    axes[1].set_title("Relation between expression and adjusted pvalue")
    axes[2].set_title("Relation between expression and cooks distance")

    # seaborn.countplot(
    #     ax=axes[0],
    #     x="Status",
    #     data=tmp
    # )
    seaborn.histplot(
        ax=axes[0],
        x="Mean Expression",
        data=tmp,
        bins=100,
        hue="Independent Filtering",
        multiple="stack"
    )
    seaborn.scatterplot(
        ax=axes[1],
        data=tmp,
        y="-Log10(P-Adjusted)",
        x="Log10(Mean Expression)",
        hue="Independent Filtering",
        style="Status"
    )
    seaborn.scatterplot(
        ax=axes[2],
        data=tmp,
        y="Maximum Cook's distance",
        x="Log10(Mean Expression)",
        hue="Independent Filtering",
        style="Status"
    )

    matplotlib.pyplot.savefig(
        png_out,
        bbox_inches="tight"
    )

    matplotlib.pyplot.close()


def pval_distribution(deseq: pandas.DataFrame,
                      comparison: str,
                      png_out: str,
                      chromosomes: Optional[List[str]] = None) -> None:
    """Plot number of differentially expressed genes and their Chromosomes"""
    print(deseq.columns.tolist())
    tmp = deseq[
        ["padj", "Status", "filterThreshold", "Chromosome"]
    ]
    tmp["padj"] = -numpy.log10(tmp["padj"])
    tmp["filterThreshold"] = [
        "Pass" if filter is True else "Do not pass"
        for filter in tmp["filterThreshold"]
    ]

    if chromosomes is not None:
        tmp = tmp[tmp["Chromosome"].isin(chromosomes)]

    tmp.sort_values(by="Chromosome", inplace=True)

    tmp.columns = [
        "-Log10(P-Adjusted)",
        "Status",
        "Independent Filtering",
        "Chromosomes"
    ]

    # Build plot
    fig, axes = matplotlib.pyplot.subplots(
        1, 2, figsize=(20,10)#, sharey=True
    )
    fig.suptitle(f"Differential gene expression {comparison}")
    axes[0].set_title("# differentially expressed genes")
    axes[1].set_title("# differentially expressed genes per chromosomes")

    seaborn.countplot(
        ax=axes[0],
        x="Status",
        data=tmp
    )

    seaborn.countplot(
        ax=axes[1],
        x="Status",
        data=tmp,
        hue="Chromosomes"
    )

    matplotlib.pyplot.savefig(
        png_out,
        bbox_inches="tight"
    )
    matplotlib.pyplot.close()


def filter_theta(theta: pandas.DataFrame,
                 metadata: pandas.DataFrame,
                 png_out: str) -> None:
    """Plot number of rejected observations due to theta filer"""
    theta.columns = ["Quantiles of filter", "# of rejections"]
    fig = matplotlib.pyplot.figure(figsize=(10,10))
    seaborn.lineplot(
        data=theta,
        x="Quantiles of filter",
        y="# of rejections"
    )
    matplotlib.pyplot.axvline(
        x=metadata["filterTheta"].tolist()[0],
        color="red",
        linestyle="--"
    )
    matplotlib.pyplot.savefig(
        png_out,
        bbox_inches="tight"
    )

    matplotlib.pyplot.close()


# Setting style and colors
seaborn.set_theme(style="whitegrid")
condition_list = list(snakemake.params["condition_dict"].values())
colors = seaborn.husl_palette(len(condition_list), s=0.45)
cond_colors = {
    str(cond): color for cond, color in zip(condition_list, list(colors))
}
sample_colors = {
    s: cond_colors[snakemake.params["condition_dict"][s]]
    for s in snakemake.params["condition_dict"].keys()
}

# Loading data
deseq, counts, mu, theta, metadata = load_all_data(
    deseq_path=snakemake.input["deseq"],
    intermediar_path=snakemake.input["intermediar"],
    dst_path=snakemake.input["dst"],
    assays_path=snakemake.input["assays"],
    filter_theta_path=snakemake.input["filter_theta"],
    metadata_path=snakemake.input["metadata"],
    gene2gene_path=snakemake.input["gene2gene"],
    padj_threshold=snakemake.params.get("padj_threshold", 0.05)
)

# Plotting
if "log_counts" in snakemake.output.keys():
    plot_counts(
        dst=counts,
        colors=sample_colors,
        png_out=snakemake.output["log_counts"]
    )

if "log_mu" in snakemake.output.keys():
    plot_counts(
        dst=mu,
        colors=sample_colors,
        png_out=snakemake.output["log_mu"],
        is_mu=True
    )

if "gene_plots" in snakemake.output.keys():
    try:
        gene_list = snakemake.wildcards["gene"]
    except AttributeError:
        gene_list = snakemake.params["gene_list"]

    plot_single_gene(
        dst=counts,
        deseq=deseq,
        mu=mu,
        condition_dict=snakemake.params["condition_dict"],
        genes=gene_list,
        png_prefix=snakemake.params["gene_plots_prefix"],
        comparison=snakemake.params["comparison"]
    )

if "pval" in snakemake.output.keys():
    pval_distribution(
        deseq=deseq,
        comparison=snakemake.params["comparison"],
        png_out=snakemake.output["pval"],
        chromosomes=snakemake.params.get("chromosomes", None)
    )

if "independent_filtering" in snakemake.output.keys():
    independent_filtering(
        deseq=deseq,
        comparison=snakemake.params["comparison"],
        png_out=snakemake.output["independent_filtering"]
    )

if "filter_theta" in snakemake.output.keys():
    filter_theta(
        theta=theta,
        metadata=metadata,
        png_out=snakemake.output["filter_theta"]
    )
