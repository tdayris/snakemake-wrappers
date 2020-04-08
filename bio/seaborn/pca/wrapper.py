#!/usr/bin/python3.8
# conding: utf-8

"""
Plot a PCA
"""

import itertools
import logging
import matplotlib
import matplotlib.pyplot
import numpy
import pandas
import seaborn
import sklearn.decomposition


from os.path import basename, dirname, commonprefix
from snakemake.utils import makedirs

logging.basicConfig(
    filename=snakemake.log[0],
    filemode="w",
    level=logging.DEBUG
)

# Load data and remove text annotations
data = pandas.read_csv(
    snakemake.input["counts"],
    sep="\t",
    header=0,
    index_col=0
)
data = data[list(data.select_dtypes(include=[numpy.number]).columns.values)]

condition_dict = snakemake.params.conditions

# Perform PCA
nbc = len(data.columns.tolist())
skpca = sklearn.decomposition.PCA(n_components=nbc)

# Prepare plots
sktransform = skpca.fit_transform(data.T)
skvar = skpca.explained_variance_ratio_
results = pandas.DataFrame(
    sktransform,
    columns=[f"PC{i}" for i in range(1, nbc+1, 1)],
    index=data.columns.tolist()
)
logging.debug(results)

seaborn.set(style="darkgrid")
output_prefix = snakemake.params.get("prefix", "pca")
legend_position = snakemake.params.get("legend_position", "upper center")
axes = snakemake.params.get("axes", range(1, 4, 1))
for ax1, ax2 in itertools.permutations(axes, 2):
    results["Conditions"] = [
        condition_dict[i] for i in results.index
    ]

    g = seaborn.FacetGrid(
        results,
        hue="Conditions",
        height=13
    )

    name_ax1 = f"PC{ax1}"
    skvar_ax1 = skvar[ax1 - 1] * 100
    name_ax2 = f"PC{ax2}"
    skvar_ax2 = skvar[ax2 - 1] * 100
    logging.info(f"Building plot: {output_prefix}_{name_ax1}_{name_ax2}.png")

    g = g.map(
        matplotlib.pyplot.scatter,
        name_ax1,
        name_ax2,
        s=50
    )

    matplotlib.pyplot.title(
        f"{name_ax1} ({skvar_ax1:.2f}) and {name_ax2} (skvar_ax2:.2f)"
    )

    if snakemake.params.get("samples_names", False) is True:
        points_coordinates = zip(
            results.index.tolist(),
            results[name_ax1],
            results[name_ax2]
        )

        for label, x, y in points_coordinates:
            matplotlib.pyplot.annotate(
                label,
                xy=(x, y),
                xytext=(-5, -5),
                textcoords="offset points",
                ha="center",
                va="top"
            )

    frame = (matplotlib.pyplot
                       .legend(loc=legend_position, frameon=True)
                       .get_frame())
    frame.set_facecolor("white")

    matplotlib.pyplot.savefig(
        f"{output_prefix}_{name_ax1}_{name_ax2}.png",
        bbox_inches="tight"
    )
logging.info("Process over")
