#!/usr/bin/python3.8
# coding: utf-8

""" Snakemake wrapper for pair grid """
import logging

import seaborn
import logging
import matplotlib.pyplot
import pandas


def read_tsv(path: str, 
             key: str = "ID", 
             values: list[str] = ["Count", "p.adjust", "GeneRatio"]
            ) -> pandas.DataFrame:
    """Return a dataframe from a single tsv file"""
    logging.debug("Loading %s", path)
    tmp = pandas.read_csv(
        path,
        sep="\t",
        header=0,
        index_col=key
    )[[values]]


def load_tsv(paths: list[str],  
             key: str = "ID",
             values: list[str] = ["Count", "p.adjust", "GeneRatio"]
            ) -> pandas.DataFrame:
    """Return a dataframe from multiple TSV files"""
    result = None
    for path in paths:
        tmp = read_tsv(path, key, values)
        try:
            result = pandas.merge(
                left=result,
                right=tmp,
                left_index=True,
                right_index=True,
            )
        except AttributeError:
            result = tmp
    logging.debug("All tables loaded: %s", str(result.head()))
    return result

sns.set_theme(style="whitegrid")

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


# Load user input depending on their number
df = None
if isinstance(snakemake.input["tsv"], str):
    df = read_tsv(snakemake.input["tsv"], key=snakemake.params.get("key"))
else:
    df = load_tsv(snakemake.input["tsv"], key=snakemake.params.get("key"))


# Make the PairGrid
g = seaborn.PairGrid(
    df.sort_values(snakemake.params["x_axis"], ascending=False),
    x_vars=df.columns,
    y_vars=[snakemake.params["y_axis"]],
    height=10,
    aspect=.25
)

# Draw a dot plot using the stripplot function
g.map(seaborn.stripplot, size=10, orient="h", jitter=False,
      palette="flare_r", linewidth=1, edgecolor="w")

# Use the same x axis limits on all columns and add better labels
g.set(xlim=(0, 25), xlabel=snakemake.params["y_axis"], ylabel="")

# Use semantically meaningful titles for the columns

for ax, title in zip(g.axes.flat, snakemake.params["subplots"]):

    # Set a different title for each axes
    ax.set(title=title)

    # Make the grid horizontal instead of vertical
    ax.xaxis.grid(False)
    ax.yaxis.grid(True)

# Save figure
seaborn.despine(left=True, bottom=True)
matplotlib.pyplot.savefig(
    snakemake.output["png"],
    bbox_inches="tight"
)