#!/usr/bin/python3.8
# conding: utf-8

"""
Plot a clustered heatmap of a set of genes
"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2020, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

import logging
import matplotlib
import matplotlib.pyplot
import numpy
import pandas
import seaborn

from os.path import basename, dirname
from snakemake.utils import makedirs

logging.basicConfig(
    filename=snakemake.log[0],
    filemode="w",
    level=logging.DEBUG
)


# Build output directory if necessary
if (outdir := basename(dirname(snakemake.output["png"]))) != "":
    makedirs(outdir)
    logging.debug(f"Directory: '{outdir}' created.")

# Load normalized counts
data = pandas.read_csv(
    snakemake.input["tsv"],
    sep="\t",
    header=0,
    index_col=4
)
logging.debug("Head of the loaded data:")
logging.debug(data.head())

# Keep samples anntoated in "conditions" parameter
conditions = snakemake.params["conditions"]
sample_list = list(conditions.keys())
condition_set = set(conditions.values())
# data.set_index("Gene_Name", inplace=True)
data = data[sample_list]

# Create custom colormap for heatmap values
cmap = seaborn.diverging_palette(
    h_neg=240,
    h_pos=10,
    as_cmap=True
)
logging.info("Color palette built")

# Create a categorical palette for samples identification
colors = seaborn.husl_palette(len(condition_set), s=0.45)
cond_colors = {
    str(cond): color for cond, color in zip(list(condition_set), list(colors))
}
sample_colors = {
    sample: cond_colors[conditions[sample]]
    for sample in sample_list
}
logging.info(f"Color palette assigned to samples: {sample_colors}")
conditions = pandas.DataFrame.from_dict(conditions, orient="index")
cond_id = snakemake.params.get("factor", "Condition")
conditions.columns = [cond_id]
logging.info(conditions.head())

# Sorry for that part, yet I could not manage to find any
# other way to perform quicker multi-level indexing
data = data.T
data = pandas.merge(
    data,
    conditions,
    left_index=True,
    right_index=True,
    how="left"
)
data = (data.reset_index()
            .set_index([cond_id, "index"])
            .T)
logging.info("Multi level indexing build")

condition_colors = (
    pandas.Series(data.columns.get_level_values(cond_id), index=data.columns)
          .map(cond_colors)
)
logging.info("Sample/Color mapped")

# Build graph
# data = data.corr()
print(len(data.columns.tolist()) / 50)
ax = seaborn.clustermap(
    data,
    cmap=cmap,
    col_colors=(
        condition_colors
        if snakemake.params.get("col_condition_color", True) is True
        else None
    ),
    method="average",
    metric="euclidean",
    # zscore=0,
    standard_scale=0,
    row_cluster=(snakemake.params.get("row_cluster", True) is True),
    col_cluster=(snakemake.params.get("col_cluster", True) is True),
    linewidths=0,
    figsize=(
        7 if len(data.columns.tolist()) <= 20 else 15,
        7 if len(data.columns.tolist()) <= 100 else 25
    ),
    robust=(snakemake.params.get("robust", False) is True)
)
logging.info("clustermap built")

# Rotate sample id to make them readable
matplotlib.pyplot.setp(
    ax.ax_heatmap.yaxis.get_majorticklabels(),
    rotation=snakemake.params.get("ylabel_rotation", 0)
)

matplotlib.pyplot.setp(
    ax.ax_heatmap.xaxis.get_majorticklabels(),
    rotation=snakemake.params.get("xlabel_rotation", 90)
)

# Save result
matplotlib.pyplot.savefig(
    snakemake.output["png"],
    bbox_inches="tight"
)
logging.info("Process over")
