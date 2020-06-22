#!/usr/bin/python3.8
# conding: utf-8

"""
This script computes relative position of the variant in its parent gene.
"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2020, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

import logging
import re

vcf_path = snakemake.input["vcf"]
tsv_path = snakemake.output["tsv"]
gene_re = re.compile(snakemake.params.get("gene_re", "ENSG[0-9]+"))

logging.basicConfig(
    filename=snakemake.log[0],
    filemode="w",
    level=logging.DEBUG
)

with open(vcf_path, "r") as vcf, open(tsv_path, "w") as tsv:
    # Write header on user's request
    if snakemake.params.get("header", False) is True:
        cols = ["Chromosome", "GeneID", "Position", "Ref", "Alt\n"]
        tsv.write('\t'.join(cols))

    for line in vcf:
        if line.startswith("#"):
            # Then this is a comment
            continue

        chomp = line[:-1].split("\t")
        genes_id = gene_re.findall(chomp[7])

        if len(genes_id) == 0:
            logging.warning(f"WARNING: no gene id found in {line}")
            continue

        for gene in genes_id:
            tsv.write(
                "\t".join([chomp[0], chomp[1], gene, chomp[2], chomp[3]])
            )
            tsv.write("\n")
