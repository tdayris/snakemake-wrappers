#!/usr/bin/python3.8
# conding: utf-8

"""
Filter a DESeq2 output file
"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2020, Thibault Dayris"
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
    snakemake.input["tsv"],
    sep="\t",
    header=0,
    index_col=0,
    dtype={
        1: str,
        2: numpy.float,
        3: numpy.float,
        4: numpy.float,
        5: numpy.float,
        6: numpy.float,
        7: numpy.float,
    }
)

print(data.head())

if (gene2gene := snakemake.input.get("gene2gene", None)) is not None:
    genetable = pandas.read_csv(
        gene2gene,
        sep="\t",
        header=0,
        index_col=None
    )

    data = pandas.merge(
        data.copy(),
        genetable,
        left_index=True,
        right_on="Gene_ID",
        how="left"
    )

    data = data[[
        "Gene_ID", "Gene_Name", "log2FoldChange", "padj",
        "Chromosome", "Start", "Stop", "Strand"
    ]]
    data.rename(columns={"Gene_ID": "index"}, inplace=True)
else:
    data.reset_index(inplace=True)
    data = data[["index", "log2FoldChange", "padj"]]

print(data.head())

padjthreshold = float(snakemake.params.get("alpha", 0.05))
fc_threshold = float(snakemake.params.get("fold_change", 0.01))
general_table_only = snakemake.params.get("general_table_only", False)


data["Cluster_FC"] = [
    get_fc_cluster(fc, fc_threshold)
    for fc in data["log2FoldChange"]
]

data["Cluster_Sig"] = [
    get_alpha_cluster(padj, padjthreshold)
    for padj in data["padj"]
]

if "fc_sig" in snakemake.output.keys():
    logging.debug("Prining the log2(FC) / Significance table")
    tmp = data[data["Cluster_Sig"] != "Non_Significative"]
    tmp.rename(
        columns={
            "index": "GeneIdentifier",
            "log2FoldChange": "stat_change",
            "Cluster_Sig": "cluster"
        },
        inplace=True
    )
    tmp.to_csv(
        snakemake.output.fc_sig,
        sep="\t",
        index=False
    )

if "fc_fc" in snakemake.output.keys():
    logging.debug("Prining the log2(FC) / FC cluster table")
    tmp = data[data["Cluster_Sig"] != "Non_Significative"]
    tmp = tmp[tmp["Cluster_FC"] != "Non_Significative"]
    tmp.rename(
        columns={
            "log2FoldChange": "stat_change",
            "Cluster_FC": "cluster",
            "index": "GeneIdentifier"
        }
    )
    tmp.to_csv(
        snakemake.output.fc_fc,
        sep="\t",
        index=False
    )

if "padj_sig" in snakemake.output.keys():
    logging.debug("Prining the adjusted P-Value / Significance table")
    tmp = data[data["Cluster_Sig"] != "Non_Significative"]
    tmp.rename(
        columns={
            "padj": "stat_change",
            "Cluster_Sig": "cluster",
            "index": "GeneIdentifier"
        },
        inplace=True
    )
    tmp.to_csv(
        snakemake.output.padj_sig,
        sep="\t",
        index=False
    )

if "padj_fc" in snakemake.output.keys():
    logging.debug("Prining the adjusted P-Value / FoldChange table")
    tmp = data[data["Cluster_FC"] != "Non_Significative"]
    tmp.rename(
        columns={
            "padj": "stat_change",
            "index": "GeneIdentifier",
            "Cluster_FC": "cluster"
        },
        inplace=True
    )
    tmp.to_csv(
        snakemake.output.padj_fc,
        sep="\t",
        index=False
    )

if "complete" in snakemake.output.keys():
    logging.debug("Prining the complete table")
    tmp.rename(
        columns={
            "log2FoldChange": "stat_change",
            "index": "GeneIdentifier",
            "Cluster_FC": "cluster",
            "Cluster_Sig": "significance",
            "padj": "Adjuster_PValue"
        },
        inplace=True
    )
    tmp.to_csv(
        snakemake.output.complete,
        sep="\t",
        index=False
    )
