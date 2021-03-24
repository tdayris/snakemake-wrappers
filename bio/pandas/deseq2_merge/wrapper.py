#!/usr/bin/python3.8
# conding: utf-8

"""
Filter a DESeq2 output file
"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2021, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

import logging
import pandas
import numpy

def get_fc_cluster(value: numpy.float,
                   threshold: numpy.float = numpy.float(0.01)) -> str:
    """
    This functon returns a class for a given log2(Fold Change)
    """
    if abs(value) < threshold:
        return "Non_Significative"
    if value > 0:
        return "Up_Regulated"
    return "Down_Regulated"


def get_alpha_cluster(value: numpy.float,
                      threshold: numpy.float = numpy.float(0.05)) -> str:
    """
    This function returns a class for a given adjusted pval
    """
    if value < threshold:
        return "Differentially_Expressed"
    return "Non_Significative"

logging.basicConfig(
    filename=snakemake.log[0],
    filemode="w",
    level=logging.DEBUG
)

data = pandas.read_csv(
    snakemake.input["wald_tsv"],
    sep="\t",
    header=0,
    index_col=0,
    dtype={
        1: str,
        2: float,
        3: float,
        4: float,
        5: float,
        6: float,
        7: float,
    }
)

logging.debug(data.head())

genetable = pandas.read_csv(
    snakemake.input.get("gene2gene"),
    sep="\t",
    header=0,
    index_col=None
)
genetable.set_index("Gene_ID", inplace=True)
logging.debug(genetable.head())

data = pandas.merge(
    data.copy(),
    genetable,
    left_index=True,
    right_index=True,
    how="left"
)
logging.debug(data.head())
data = data[[
    "Gene_Name", "log2FoldChange", "padj"
]]

logging.debug(data.head())

padjthreshold = float(snakemake.params.get("alpha", 0.05))
fc_threshold = float(snakemake.params.get("fold_change", 0.001))


data["Cluster_FC"] = [
    get_fc_cluster(fc, fc_threshold)
    for fc in data["log2FoldChange"]
]

data["Cluster_Sig"] = [
    get_alpha_cluster(padj, padjthreshold)
    for padj in data["padj"]
]

counts = pandas.read_csv(
    snakemake.input["dst_tsv"],
    sep="\t",
    header=0,
    index_col=0
)

data = pandas.merge(
    data.copy(),
    counts,
    left_index=True,
    right_index=True
)
logging.debug(data.head())
try:
    data.set_index("Gene_Name", inplace=True)
except IndexError:
    logging.debug(
        "Duplicate gene name in Gene_Name column. They won't be used as genes "
        "identifiers since it would lead to errors."
    )
logging.debug(data.head())

if "filtered_counts" in snakemake.output.keys():
    logging.debug("Saving filtered counts to TSV")
    deseq_cols = {
        "log2FoldChange", "padj", "Cluster_Sig",
        "Cluster_FC", "index", "Gene_Name"
    }
    samples_id = set(data.columns.tolist()) - deseq_cols
    tmp = data[data["Cluster_Sig"] == "Differentially_Expressed"]
    tmp = tmp[tmp["Cluster_FC"] != "Non_Significative"]
    tmp = tmp[list(samples_id)]
    tmp.to_csv(snakemake.output["filtered_counts"], sep="\t")
    del tmp

if "filtered_deseq2" in snakemake.output.keys():
    logging.debug("Savig filtered deseq2 to TSV")
    deseq_cols = [
        "log2FoldChange", "padj"
    ]
    tmp = data[data["Cluster_Sig"] == "Differentially_Expressed"]
    tmp = tmp[tmp["Cluster_FC"] != "Non_Significative"]
    tmp = tmp[deseq_cols]
    tmp.to_csv(snakemake.output["filtered_deseq2"], sep="\t")
    del tmp

if "merged_table" in snakemake.output.keys():
    logging.debug("Saving complete table to TSV")
    data.to_csv(snakemake.output["merged_table"], sep="\t")
    logging.debug(data.head())
