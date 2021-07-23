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

conditions = pandas.DataFrame.from_dict(
    snakemake.params["conditions"],
    orient="index"
)
cond_id = snakemake.params.get("factor", "Condition")
conditions.columns = [cond_id]
logging.debug("Head of the conditions dataframe:")
logging.debug(conditions.head())

# Load normalized counts
data = pandas.read_csv(
    snakemake.input["counts"],
    sep="\t",
    header=0,
    index_col=0
)
logging.debug("Head of the loaded data:")
logging.debug(conditions.head())

# Remove possible text annotations and validate
data = data[list(data.select_dtypes(include=[numpy.number]).columns.values)]

# Create custom colormap for heatmap values
cmap = seaborn.diverging_palette(
    h_neg=240,
    h_pos=10,
    as_cmap=True
)
logging.info("Color palette built")

# Create a categorical palette for samples identification
cond_set = set(conditions[cond_id])
colors = seaborn.husl_palette(len(cond_set), s=0.45)
cond_colors = {
    str(cond): color for cond, color in zip(list(cond_set), list(colors))
}
sample_colors = {
    sample: cond_colors[cond]
    for sample, cond in zip(conditions.index, conditions[cond_id])
}
logging.info("Color palette assigned to samples")

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

condition_colors = (pandas.Series(data.columns.get_level_values(cond_id),
                                  index=data.columns)
                          .map(cond_colors))
logging.info("Sample/Color mapped")

# Build graph
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
    zscore=0,
    row_cluster=(snakemake.params.get("row_cluster", True) is True),
    col_cluster=(snakemake.params.get("col_cluster", True) is True),
    linewidths=0.5,
    figsize=(
        min(15, int(len(data.columns.tolist())/50) * 10),
        min(15, int(len(data.index.tolist())/50) * 10)
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
