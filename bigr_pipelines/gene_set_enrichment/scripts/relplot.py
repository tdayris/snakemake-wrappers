#!/usr/bin/env python3
# coding: utf-8

"""Plot pairwise stirplot from clusterprofiler results"""

import logging
import matplotlib.pyplot
import numpy
import pandas
import seaborn


def compute_ratio(ratio: str) -> float:
    """Based on string ratio, compute its value"""
    num, den = map(float, ratio.split("/"))
    return num / den


def min_log10(x: float) -> float:
    """Compute -log10(x)"""
    return - numpy.log10(x)


def read_clusterprofiler(path: str,
                         id_col: str,
                         count_col: str,
                         padj_col: str,
                         ratio_col: str,
                         name: str
                        ) -> pandas.DataFrame:
    """Load ClusterProfiler in memory"""
    logging.debug("Loading %s in memory", path)
    tmp = pandas.read_csv(
        path,
        header=0,
        sep=" ",
        index_col=None,
    )
    tmp = tmp[[id_col, count_col, padj_col, ratio_col]]
    tmp["Comparisons"] = [name for _ in tmp.index]
    tmp[f"-log10({padj_col})"] = tmp[padj_col].apply(min_log10)
    tmp["GeneRatio"] = tmp[ratio_col].apply(compute_ratio)
    
    logging.debug(tmp.head())
    return tmp


# Set logging behaviour
try:
    # Case user provided logging file
    logging.basicConfig(
        filename=snakemake.log[0],
        filemode="w",
        level=logging.DEBUG
    )
except IndexError:
    # Case user did not provide any logging file
    logging.basicConfig(filemode="w", level=logging.DEBUG)
logging.getLogger('matplotlib').setLevel(logging.WARNING)

# Tell me how to read cluster profiler
id_col = snakemake.params.get("id_col", "ID")
count_col = snakemake.params.get("counts", "Count")
padj_col = snakemake.params.get("padj", "p.adjust")
ratio_col = snakemake.params.get("gene_ratio", "GeneRatio")

# Tell me the name of the comparisons
input_iterator = zip(snakemake.input["tsv"], snakemake.params["comparisons"])

frames = []
for path, comparison in input_iterator:
    frames.append(
        read_clusterprofiler(
            path=path,
            id_col=id_col,
            count_col=count_col,
            padj_col=padj_col,
            ratio_col=ratio_col,
            name=comparison
        )
    )

logging.debug("Building complete dataframe")
df = pandas.concat(frames, join="inner")
logging.debug(df.columns.tolist())
logging.info(df.head())

logging.debug("Building graph")
seaborn.set_theme(style="whitegrid")

g = seaborn.relplot(
    x=id_col,
    y="GeneRatio",
    hue="Comparisons",
    size=f"-log10({padj_col})",
    legend="brief",
    kind="scatter",
    data=df
)

logging.debug("Saving figure")
matplotlib.pyplot.savefig(
    snakemake.output[0],
    bbox_inches="tight"
)