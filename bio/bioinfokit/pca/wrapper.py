"""Snakemake wrapper for bioinfokit PCA"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2020"
__email__ = "Thibault Dayris"
__license__ = "MIT"

import bioinfokit.visuz
import csv
import logging
import matplotlib.pyplot
import numpy
import os
import pandas
import seaborn

from snakemake.shell import shell
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from tempfile import TemporaryDirectory


logging.basicConfig(
    filename=snakemake.log[0],
    filemode="w",
    level=logging.DEBUG
)
logging.getLogger('matplotlib.font_manager').disabled = True

# Detect delimiter with python
with open(snakemake.input[0], "r") as table:
    dialect = csv.Sniffer().sniff(table.readline())

# Loading dataset
counts = pandas.read_csv(
    snakemake.input[0],
    sep=dialect.delimiter,
    **snakemake.params.get("read_csv", {})
)

if snakemake.params.get("standardize", False) is True:
    logging.debug("Standardizing data")
    tmp = StandardScaler().fit_transform(counts)
    counts = pandas.DataFrame(tmp, columns=counts.columns)
    del tmp

logging.debug(counts.head())

pca = PCA().fit(counts)
pca_scores = PCA().fit_transform(counts)
loadings = pca.components_
num_pc = pca.n_features_
variance_ratio = pca.explained_variance_ratio_

pc_list = ["PC"+str(i) for i in list(range(1, num_pc+1))]
loadings_df = pandas.DataFrame.from_dict(dict(zip(pc_list, loadings)))
loadings_df['variable'] = counts.columns.values
loadings_df = loadings_df.set_index('variable')
logging.debug("PCA loadings:")
logging.debug(loadings_df)

if "loadings_correlation_heatmap" in snakemake.output.keys():
    logging.debug("Saving correlation heatmap")
    seaborn.clustermap(
        loadings_df, annot=True, cmap='Spectral',
        row_cluster=True, col_cluster=True
    )
    matplotlib.pyplot.savefig(
        snakemake.output["loadings_correlation_heatmap"],
        bbox_inches="tight"
    )
    matplotlib.pyplot.clf()


if "pca_scree" in snakemake.output.keys():
    with TemporaryDirectory() as tempdir:
        logging.debug("Saving PCA Screen plot")
        figtype = snakemake.output["pca_scree"].split(".")[-1]
        bioinfokit.visuz.cluster.screeplot(
            obj=[pc_list, variance_ratio],
            figtype=figtype,
            show=False
        )
        os.rename("screeplot.png", snakemake.output["pca_scree"])


if "pca_2d" in snakemake.output.keys():
    logging.debug("Saving 2D PCA")
    figtype = snakemake.output["pca_2d"].split(".")[-1]
    bioinfokit.visuz.cluster.pcaplot(
        x=loadings[0],
        y=loadings[1],
        figtype=figtype,
        labels=counts.columns.values,
        var1=round(variance_ratio[0]*100, 2),
        var2=round(variance_ratio[1]*100, 2),
        show=False
    )
    os.rename("pcaplot_2d.png", snakemake.output["pca_2d"])


if "pca_3d" in snakemake.output.keys():
    logging.debug("Saving 3D PCA")
    figtype = snakemake.output["pca_3d"].split(".")[-1]
    bioinfokit.visuz.cluster.pcaplot(
        x=loadings[0],
        y=loadings[1],
        z=loadings[2],
        figtype=figtype,
        labels=counts.columns.values,
        var1=round(variance_ratio[0]*100, 2),
        var2=round(variance_ratio[1]*100, 2),
        var3=round(variance_ratio[2]*100, 2),
        show=False
    )
    os.rename("pcaplot_3d.png", snakemake.output["pca_3d"])


if "biplot_2d" in snakemake.output.keys():
    logging.debug("Saving 2D biplot")
    figtype = snakemake.output["biplot_2d"].split(".")[-1]
    bioinfokit.visuz.cluster.biplot(
        cscore=pca_scores,
        loadings=loadings,
        figtype=figtype,
        labels=counts.columns.values,
        var1=round(variance_ratio[0]*100, 2),
        var2=round(variance_ratio[1]*100, 2),
        show=False
    )
    os.rename("biplot_2d.png", snakemake.output["biplot_2d"])

if "biplot_3d" in snakemake.output.keys():
    logging.debug("Saving 3D biplot")
    figtype = snakemake.output["biplot_3d"].split(".")[-1]
    bioinfokit.visuz.cluster.biplot(
        cscore=pca_scores,
        loadings=loadings,
        figtype=figtype,
        labels=counts.columns.values,
        var1=round(variance_ratio[0]*100, 2),
        var2=round(variance_ratio[1]*100, 2),
        var3=round(variance_ratio[2]*100, 2),
        show=False
    )
    os.rename("biplot_3d.png", snakemake.output["biplot_3d"])
