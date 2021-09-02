#!/usr/bin/python3.8
# conding: utf-8

"""
Annotate MAF file with Cancer sensus genes
"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2021, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"


import csv
import pandas
import numpy

def read_cancer_sensus(path: str) -> pandas.DataFrame:
    """
    Load official cancer gene sensus TSV/CSV in memory
    """
    # Detect delimiter with python
    with open(path, "r") as table:
        dialect = csv.Sniffer().sniff(table.readline())
        logging.debug(f"Detected dialect: {dialect.delimiter}")

    # Loading dataset
    return pandas.read_csv(
        path,
        sep=dialect.delimiter,
        index_col=None,
        header=0,
        dtype={
            "Gene Symbol": str,
            "Name": str,
            "Entrez GeneId": int,
            "Genome Location": str,
            "Tier": int,
            "Hallmark": str,
            "Chr Band": numpy.float,
            "Somatic": "category",
            "Germline": "category",
            "Tumour Types(Somatic)": str,
            "Tumour Types(Germline)": str,
            "Cancer Syndrome": str,
            "Tissue Type": "category",
            "Molecular Genetics": "category",
            "Role in Cancer": "category",
            "Mutation Types": "category",
            "Translocation Partner": str,
            "Other Germline Mut": str,
            "Other Syndrome": str,
            "Synonyms": str,

        }
    )


def read_maf(path: str) -> pandas.DataFrame:
    return pandas.read_csv(
        path,
        sep=dialect.delimiter,
        index_col=None,
        header=0
    )


# Loading input files
sanger = read_cancer_sensus(snakemake.input["sanger"])
logging.info("Sanger cancer gene sensus dataset loaded")
logging.debug(sanger.head())

maf = read_maf(snakemake.input["maf"])
logging.info("MAF file loaded")
logging.debug(maf.head())

# Merging on gene name
maf = pandas.merge(
    left=maf,
    right=sanger,
    how="left",
    left_on=snakemake.params.get("left_key", ""),
    right_on=snakemake.params.get("right_key", ""),
    suffixes=("", "_SangerCGC")
)
logging.info("Merge performed")
logging.debug(maf.head())

maf.to_csv(
    snakemake.output["maf"],
    sep="\t",
    index=True,
    header=True
)
logging.info("Process over")
