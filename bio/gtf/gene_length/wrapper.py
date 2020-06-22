#!/usr/bin/python3.8
# conding: utf-8

"""
This script extracts a tsv from a GTF. This tsv contains gene_id
and gene_name and gene length. Optionaly, headers and positions can
also be written.
"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2020, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

gtf_in = snakemake.input["gtf"]
tsv_out = snakemake.output["tsv"]

with open(gtf_in, "r") as gtf, open(tsv_out, "w") as tsv:
    # Write header on user request
    if snakemake.params.get("header", False) is True:
        cols = (
            ["Gene_ID", "Gene_Name", "Length", "Chromosome",
             "Start", "Stop", "Strand"]
            if snakemake.params.get("positions", False) is True else
            ["Gene_ID", "Gene_Name", "Length"]
        )
        tsv.write('\t'.join(cols))
        tsv.write("\n")

    for line in gtf:
        if line.startswith("#"):
            # Then is it a comment
            continue

        # If the line is not a transcript, then it will no contain
        # the required information
        chomp = line[:-1].split("\t")
        if chomp[2] != "gene":
            continue

        start = chomp[3]
        stop = chomp[4]
        chrom = chomp[0]
        strand = chomp[6]

        # We are interested in the last column
        chomp = {
            attr.split('"')[0].strip(): attr.split('"')[1].strip()
            for attr in chomp[8].split(";")
            if attr != '' and '"' in attr
        }

        if snakemake.params.get("gencode", False) is True:
            chomp["gene_id"] = chomp["gene_id"].split(".")[0]


        # Some genes have an ID but no name ...
        try:
            result = [
                chomp["gene_id"],
                chomp["gene_name"],
                str(abs(int(stop) - int(start)))
            ]
        except KeyError:
            result = [
                chomp["gene_id"],
                chomp["gene_id"],
                str(abs(int(stop) - int(start)))
            ]

        if snakemake.params.get("positions", False) is True:
            tsv.write("\t".join(result + [chrom, start, stop, strand]) + "\n")
        else:
            tsv.write("\t".join(result) + "\n")
