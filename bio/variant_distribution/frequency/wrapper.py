#!/usr/bin/python3.8
# conding: utf-8


"""
This script takes two tables built by previous wrappers and builds a
variant frequency table for future plots.
"""


import logging
import pandas
import numpy

from typing import List


def compute_count(data_frame: pandas.DataFrame,
                  group_cols: List[str],
                  count: str) -> pandas.DataFrame:
    """
    Compute value frequency in a DataFrame, based on a list of columns
    that are to be grouped together
    """
    return variants.groupby(by=group_cols)[count].transform("count")


logging.basicConfig(
    filename=snakemake.log[0],
    filemode="w",
    level=logging.DEBUG
)


# Loading files and parameters
genes_lengths = pandas.read_csv(
    snakemake.input["gene_length"],
    sep="\t",
    header=0,
    index_col=None,
    dtype={
        0: str,  # Gene ID
        1: str,  # Gene Name
        2: numpy.int,  # Length
        3: str,  # Chromosome
        4: numpy.int,  # Start
        5: numpy.int  # Stop
    }
)

logging.info("Loading datasets")
variant_positions = pandas.read_csv(
    snakemake.input["variant_positions"],
    sep="\t",
    header=0,
    index_col=None,
    dtype={
        0: str,  # Chromosome
        1: str,  # Gene ID
        2: numpy.int,  # Position
        3: str,  # Reference
        4: str  # Alternative
    }
)

round_digits = snakemake.params.get("round_digits", 2)
logging.info("Datasets loaded")


variants = pandas.merge(
    genes_lengths,
    variant_positions,
    on=["GeneID", "Chromosome"]
)
logging.info("Variant table merged with annotation:")
logging.info(variants.head())

# Filtering out annotation errors
logging.debug(variants.shape)
variants = variants[variants["Position"] > variants["Start"]]
variants = variants[variants["Position"] < variants["Stop"]]
logging.debug(variants.shape)
logging.info("Variants which position is out of their gene's bounds are filtered out.")

# Building variants positions relatively to their gene start
variant_iterator = zip(
    variants.Start,
    variants.Position,
    variants.Length
)

variants["Relative"] = [
    (position - start) / length
    for start, position, length in variant_iterator
]
logging.info("Relative position of variants computed:")
logging.debug(variants.head())

variants["RoundedRelative"] = variants.Relative.round(
    decimals=round_digits
)
logging.info("Frequency of the variants relative position computed:")
logging.debug(variants.head())

variants["GeneCount"] = compute_count(
    variants.copy(),
    ["Chromosome", "GeneID", "GeneName"],
    "GeneID"
)


variants["VariantCount"] = compute_count(
    variants.copy(),
    ["Chromosome", "GeneID", "GeneName", "RoundedRelative"],
    "RoundedRelative"
)
variants["VariantFrequency"] = [
    count / total
    for count, total in zip(variants.VariantCount, variants.GeneCount)
]
logging.info("Per-gene variant-frequency computed:")
logging.debug(variants.head())


variants["AltCount"] = compute_count(
    variants.copy(),
    ["Chromosome", "GeneID", "GeneName", "Alt", "RoundedRelative"],
    "RoundedRelative"
)
variants["AltFrequency"] = [
    count / total
    for count, total in zip(variants.AltCount, variants.GeneCount)
]
logging.info("Per-gene, allele-wise variant-frequency computed:")
logging.debug(variants.head())

variants["GeneCount"] = compute_count(
    variants.copy(),
    ["Chromosome", "GeneID", "GeneName"],
    "GeneID"
)



# del variants["RoundedRelative"]
# del variants["AltCount"]
# del variants["VariantCount"]
logging.debug(variants.head())
variants.to_csv(
    snakemake.output["tsv"],
    sep="\t",
    index=False
)
logging.info("Process over")
