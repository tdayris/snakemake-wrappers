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
    #filename=snakemake.log[0],
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

logging.debug(data.head())
base_cols = list(
    set(["log2FoldChange", "padj", "pvalue"] + snakemake.params.get("keep_cols", []))
)

if (gene2gene := snakemake.input.get("gene2gene", None)) is not None:
    logging.info("Loading gene information")
    genetable = pandas.read_csv(
        gene2gene,
        sep="\t",
        header=0,
        index_col=None
    )
    logging.debug(genetable.head())

    data = pandas.merge(
        data.copy(),
        genetable,
        left_index=True,
        right_on="Gene_ID",
        how="left"
    )
    
    base_cols += [
        "Gene_Name",
        "Chromosome", "Strand"
    ]
    logging.debug(data.head())

    data = data[["Gene_ID"] + base_cols]
    data.rename(columns={"Gene_ID": "index"}, inplace=True)
else:
    data.reset_index(inplace=True)
    data = data[["index"] + base_cols]

if (counts := snakemake.input.get("dst", None)) is not None:
    logging.info("Loading gene counts")
    counts_table = pandas.read_csv(
        counts,
        sep="\t",
        header=0,
        index_col=0
    )
    base_cols += counts_table.columns.to_list()
    logging.debug(counts_table.head())
    
    data = pandas.merge(
        data.copy(),
        counts_table,
        left_on="index",
        right_index=True,
        how="left"
    )
    
    ref_samples = snakemake.params.get("ref_samples", None)
    test_samples = snakemake.params.get("test_sample", None)
    if all(isinstance(i, list) for i in [ref_samples, test_samples]):
        base_cols += ["ratio"]
        data["ref_mean"] = data[ref_samples].mean(axis=1)
        data["test_mean"] = data[test_samples].mean(axis=1)
        data["ratio"] = data["test_mean"] / data["ref_mean"]
    
    data[["index"] + base_cols]
    

logging.debug(data.head())

if "counts" in snakemake.input.keys():
    logging.info("Loading counts")
    counts = pandas.read_csv(
        snakemake.input["counts"],
        sep="\t",
        header=0,
        index_col=0
    )
    data = pandas.merge(
        data,
        counts,
        how="left",
        left_index=True,
        right_index=True
    )

logging.debug(data.head())

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
    logging.info("Prining the log2(FC) / Significance table")
    tmp = data.copy()
    tmp = tmp[tmp["Cluster_Sig"] != "Non_Significative"]
    tmp.dropna(
        inplace=True,
        subset=["log2FoldChange", "padj"]
    )
    tmp.rename(
        columns={
            "index": "GeneIdentifier",
            "log2FoldChange": "stat_change",
            "Cluster_Sig": "cluster"
        },
        inplace=True
    )
    logging.debug(tmp.head())
    tmp.sort_values(
        by="stat_change",
        ascending=False,
        inplace=True,
        na_position="last",
    )
    tmp.to_csv(
        snakemake.output.fc_sig,
        sep="\t",
        index=False
    )
    del tmp

if "fc_fc" in snakemake.output.keys():
    logging.info("Prining the log2(FC) / FC cluster table")
    tmp = data.copy()
    tmp = tmp[tmp["Cluster_Sig"] != "Non_Significative"]
    tmp = tmp[tmp["Cluster_FC"] != "Non_Significative"]
    tmp.dropna(
        inplace=True,
        subset=["log2FoldChange", "padj"]
    )
    tmp.rename(
        columns={
            "index": "GeneIdentifier",
            "log2FoldChange": "stat_change",
            "Cluster_FC": "cluster"
        },
        inplace=True
    )
    logging.debug(tmp.head())
    tmp.sort_values(
        by="stat_change",
        ascending=False,
        inplace=True,
        na_position="last",
    )
    tmp.to_csv(
        snakemake.output.fc_fc,
        sep="\t",
        index=False
    )
    del tmp

if "padj_sig" in snakemake.output.keys():
    logging.info("Prining the adjusted P-Value / Significance table")
    tmp = data.copy()
    tmp = tmp[tmp["Cluster_Sig"] != "Non_Significative"]
    tmp.dropna(
        inplace=True,
        subset=["log2FoldChange", "padj"]
    )
    tmp.rename(
        columns={
            "padj": "stat_change",
            "Cluster_Sig": "cluster",
            "index": "GeneIdentifier"
        },
        inplace=True
    )
    logging.debug(tmp.head())
    tmp.sort_values(
        by="stat_change",
        ascending=True,
        inplace=True,
        na_position="last",
    )
    tmp.to_csv(
        snakemake.output.padj_sig,
        sep="\t",
        index=False
    )
    del tmp

if "padj_fc" in snakemake.output.keys():
    logging.info("Prining the adjusted P-Value / FoldChange table")
    tmp = data.copy()
    tmp = tmp[tmp["Cluster_Sig"] != "Non_Significative"]
    tmp.dropna(
        inplace=True,
        subset=["log2FoldChange", "padj"]
    )
    tmp.rename(
        columns={
            "padj": "stat_change",
            "index": "GeneIdentifier",
            "Cluster_FC": "cluster"
        },
        inplace=True
    )
    logging.debug(tmp.head())
    tmp.sort_values(
        by="stat_change",
        ascending=True,
        inplace=True,
        na_position="last",
    )
    tmp.to_csv(
        snakemake.output.padj_fc,
        sep="\t",
        index=False
    )
    del tmp

if "complete" in snakemake.output.keys():
    logging.info("Prining the complete table")
    tmp = data.copy()
    tmp.rename(
        columns={
            "log2FoldChange": "stat_change",
            "index": "GeneIdentifier",
            "Cluster_FC": "cluster",
            "Cluster_Sig": "significance",
            "padj": "Adjusted_PValue"
        },
        inplace=True
    )
    logging.debug(tmp.head())
    tmp.sort_values(
        by="Adjusted_PValue",
        ascending=True,
        inplace=True,
        na_position="last",
    )
    tmp.to_csv(
        snakemake.output.complete,
        sep="\t",
        index=False
    )
    del tmp