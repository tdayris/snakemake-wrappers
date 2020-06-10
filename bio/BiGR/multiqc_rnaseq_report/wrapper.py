#!/usr/bin/python3.8
# -*- coding: utf-8 -*-

"""
From the given files, produce specific MultiQC yaml files that witll ship
content info MultiQC automatically
"""

import yaml

from pathlib import Path
from snakemake.shell import shell
from typing import Any, Dict


def link_data(
        section_id: str,
        section_name: str,
        description: str,
        abs_in: str,
        abs_out: str
    ) -> Dict[str, str]:
    """
    Perform the linking and builds the config dict
    """
    log = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)
    shell("ln --force -s {abs_in} {abs_out} {log} ")
    return {
        section_id: {
            "section_name": section_name,
            "description": description
        }
    }


def link_clustermap_sample(png_path: Path, output_dir: Path) -> Dict[str, Any]:
    """
    Link the provided png_path so this would be included within MultiQC
    automatically
    """
    return link_data(
        "clustermap_sample",
        "Sample clustered heatmap",
        "Per sample clustered heatmap. This plot has been build based on "
        "normalized counts from DESeq2.<br><br>Colors run from red "
        "(identical) to blue (different). Samples are named using the "
        "following scheme: condition - sample_id",
        str(png_path.absolute()),
        f"{output_dir}/clustermap_sample_mqc.png",
    )


def link_pca_axes_correlation(
        png_path: Path,
        output_dir: Path
    ) -> Dict[str, Any]:
    """
    Link the provided png_path so this would be included within MultiQC
    automatically
    """
    return link_data(
        "pca_axes_correlation",
        "PCA Axe correlations",
        "This histogram shows which factor has a higher  correlation with "
        "the first axe of the PCA. This may help to understand why a DGE "
        "returns odd results, or which factor masks another.<br><br>"
        "This does not foresee anything for the upcomming differential "
        "analysis.",
        str(png_path.absolute()),
        f"{output_dir}/pca_axes_correlation_mqc.png"
    )


def link_pairwise_scatterplot(
        png_path: Path,
        output_dir: Path
    ) -> Dict[str, Any]:
    """
    Link the provided png_path so this would be included within MultiQC
    automatically
    """
    return link_data(
        "pairwise_scatterplot",
        "Pairwise Scatterplot",
        "This is a pairwise scatterplot. It contains: (1) in diagonal: the "
        "sample names, (2) in upper half the pearson correlation between "
        "samples, (3) in lower half, the scatterplot of each genic feature "
        "against each other. <br><br>A value of 0.95 and above shows high "
        "similarities like ones observed in cell lines. A value between "
        "0.95 and 0.80 shows average similarities. A lesser value shows "
        "noticable dissimilarities.",
        str(png_path.absolute()),
        f"{output_dir}/pairwise_scatterplot_mqc.png"
    )


def link_pca_plot(png_path: Path, output_dir: Path) -> Dict[str, Any]:
    """
    Link the provided png_path so this would be included within MultiQC
    automatically
    """
    return link_data(
        "pca_plot",
        "Principal Component Analysis",
        "Principal Component Analysis (PCA) is a common method to display "
        "samples similarities and divergences.<br><br> Each point represents "
        "a sample. Two points that are close from each others, are similar to "
        "each others. Two points, that are apart from each others, are "
        "different.",
        str(png_path.absolute()),
        f"{output_dir}/pca_plot_mqc.png"
    )


def write_yaml(output_yaml: Path, data: Dict[str, Any]) -> None:
    """
    Save given dictionnary as Yaml-formatted text file
    """
    with output_yaml.open("w") as outyaml:
        yaml.dump(data, outyaml, default_flow_style=False)


mqc_config = {"custom_data": {}}
mqc_config.update(snakemake.params)
config_outpath = Path(snakemake.output["multiqc_config"])

if "clustermap_sample" in snakemake.input.keys():
    mqc_config["custom_data"].update(
        link_clustermap_sample(
            Path(snakemake.input["clustermap_sample"]),
            config_outpath.absolute().parent
        )
    )

if "pairwise_scatterplot" in snakemake.input.keys():
    mqc_config["custom_data"].update(
        link_pairwise_scatterplot(
            Path(snakemake.input["pairwise_scatterplot"]),
            config_outpath.absolute().parent
        )
    )

if "pca_plot" in snakemake.input.keys():
    mqc_config["custom_data"].update(
        link_pca_plot(
            Path(snakemake.input["pca_plot"]),
            config_outpath.absolute().parent
        )
    )

if "pca_axes_correlation" in snakemake.input.keys():
    mqc_config["custom_data"].update(
        link_pca_axes_correlation(
            Path(snakemake.input["pca_axes_correlation"]),
            config_outpath.absolute().parent
        )
    )

write_yaml(
    output_yaml=Path(snakemake.output["multiqc_config"]),
    data=mqc_config
)
