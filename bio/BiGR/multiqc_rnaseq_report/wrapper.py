#!/usr/bin/python3.8
# -*- coding: utf-8 -*-

"""
From the given files, produce specific MultiQC yaml files that witll ship
content info MultiQC automatically
"""

import yaml

from shutil import copyfile
from pathlib import Path
from snakemake.shell import shell
from typing import Any, Dict
import logging

# Snakemake loggings
log = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)

# Python loggings
logging.basicConfig(
    filename=snakemake.log[0],
    filemode="w",
    level=logging.DEBUG
)


section_dict = {
    "clustermap_sample": {
        "section": "clustermap_sample",
        "title": "Sample clustered heatmap",
        "description": (
            "Per sample clustered heatmap. This plot has been build based on "
            "normalized counts from DESeq2.<br><br>Colors run from red "
            "(identical) to blue (different). Samples are named using the "
            "following scheme: condition - sample_id"
        ),
        "plot_name": "clustermap_sample_mqc.png"
    },
    "clustermap_genes": {
        "section": "clustermap_genes",
        "title": "Heatmap of samples clustered over genes",
        "description": (
            "Per genes clustered heatmap. This plot has been build based on "
            "filtered normalized counts from DESeq2.<br><br>Colors run from red"
            " (identical) to blue (different). Samples are named using the "
            "following scheme: condition - sample_id "
        ),
        "plot_name": "clustermap_gene_mqc.png"
    },
    "pca_axes_correlation": {
        "section": "pca_axes_correlation",
        "title": "PCA Axe correlations",
        "description": (
            "This histogram shows which factor has a higher  correlation with "
            "the first axe of the PCA. This may help to understand why a DGE "
            "returns odd results, or which factor masks another.<br><br>"
            "This does not foresee anything for the upcomming differential "
            "analysis."
        ),
        "plot_name": "pca_axes_correlation_mqc.png"
    },
    "pairwise_scatterplot": {
        "section": "pairwise_scatterplot",
        "title": "Pairwise Scatterplot",
        "description": (
            "This is a pairwise scatterplot. It contains: (1) in diagonal: the "
            "sample names, (2) in upper half the pearson correlation between "
            "samples, (3) in lower half, the scatterplot of each genic feature "
            "against each other. <br><br>A value of 0.95 and above shows high "
            "similarities like ones observed in cell lines. A value between "
            "0.95 and 0.80 shows average similarities. A lesser value shows "
            "noticable dissimilarities."
        ),
        "plot_name": "pairwise_scatterplot_mqc.png"
    },
    "pca_plot": {
        "section": "pca_plot",
        "title": "Principal Component Analysis",
        "description": (
            "Principal Component Analysis (PCA) is a common method to display "
            "samples similarities and divergences.<br><br> Each point represents "
            "a sample. Two points that are close from each others, are similar to "
            "each others. Two points, that are apart from each others, are "
            "different."
        ),
        "plot_name": "pca_plot_mqc.png"
    },
    "volcanoplot": {
        "section": "volcanoplot",
        "title": "Volcano Plot",
        "description": (
            "The volcanoplot shows the repartition of genes across two "
            "dimensions: their log2(Fold Change) and the -log10(Adjusted PValue)"
            "<br><br>Each point is a gene. The higher the point, the greater "
            "is the confidence we have in the fact that it is differentially "
            "expressed according to our factor of interest. The more apart a "
            "point is from the center of the graph, the more differentially "
            "expressed it is."
        ),
        "plot_name": "volcanoplot_mqc.png"
    },
    "distro_expr": {
        "section": "distro_expr",
        "title": "Expression Distribution",
        "description": (
            "The distribution of the expression of the genes is a common quality "
            "control performed over a differential gene analysis. <br><br> We "
            "expect all samples to have similar distribution, otherwise, it would"
            " sign a possible error before the end of the normalization process."
        ),
        "plot_name": "distro_expr_mqc.png"
    },
    "ma_plot": {
        "section": "ma_plot",
        "title": "MA-plot",
        "description": (
            "The MAplot compares M (log(ratio)) to A (mean average).<br><br> "
            "This helps to highlight possible normalisation issues. We expect no"
            " direct link between these two values; the graph should remain "
            "centered around zero"
        ),
        "plot_name": "ma_plot_mqc.png"
    },
    "consensus_cluster_plus": {
        "section": "consensus_cluster_plus",
        "title": "Consensus Clustering",
        "description": (
            "This graph has been done on DESeq2 normalized counts..<br><br> "
            "This helps to highlight possible confounding effects. We expect our sample to be well separated over the factor of interest"
        ),
        "plot_name": "consensus_cluster_plus.png"
    }
}


def write_yaml(output_yaml: Path, data: Dict[str, Any]) -> None:
    """
    Save given dictionnary as Yaml-formatted text file
    """
    with output_yaml.open("w") as outyaml:
        yaml.dump(data, outyaml, default_flow_style=False)


mqc_config = {"custom_data": {}, "sp": {}, "ignore_images": False}
mqc_config.update(snakemake.params)
config_outpath = Path(snakemake.output["multiqc_config"])

for section, plot in snakemake.input.items():
    output_dir = config_outpath.absolute().parent
    output_plot = f"{config_outpath.absolute().parent}/{section}_mqc.png"
    logging.info(
        "Working on %s (%s). Saving it as %s", section, plot, output_plot
    )
    shell("ln --relative --symbolic --force {plot} {output_plot} {log}")
    mqc_config["custom_data"][section] = {
        "section_name": section_dict[section]["title"],
        "description": section_dict[section]['description']
    }
    mqc_config["sp"][section] = {"fn": output_plot}

write_yaml(
    output_yaml=Path(snakemake.output["multiqc_config"]),
    data=mqc_config
)
