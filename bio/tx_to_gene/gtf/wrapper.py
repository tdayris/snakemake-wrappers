#!/usr/bin/python3.8
# conding: utf-8

"""
This script extracts a tsv from a GTF. This tsv contains gene_id,
transcript_id, gene_name. Optionaly, headers and positions can
also be given.
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
            ["Gene_ID", "Transcript_ID", "Gene_Name",
             "Chromosome", "Start", "Stop", "Strand"]
            if snakemake.params.get("positions", False) is True else
            ["Gene_ID", "Transcript_ID", "Gene_Name"]
        )
        tsv.write('\t'.join(cols)})
        tsv.write("\n")

    for line in gtf:
        if line.startswith("#"):
            # Then is it a comment
            continue

        # If the line is not a transcript, then it will no contain
        # the required information
        chomp = line[:-1].split("\t")
        if chomp[2] != "transcript":
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

        if snakemake.params.get("gencore", False) is True:
            chomp["gene_id"] = chomp["gene_id"].split(".")[0]
            chomp["transcript_id"] = chomp["transcript_id"].split(".")[0]

        # Some genes have an ID but no name ...
        try:
            result = [
                chomp["gene_id"],
                chomp["transcript_id"],
                chomp["gene_name"]
            ]
        except KeyError:
            result = [
                chomp["gene_id"],
                chomp["transcript_id"],
                chomp["gene_id"]
            ]

        if snakemake.params.get("positions", False) is True:
            tsv.write("\t".join(result + [chrom, start, stop, strand]) + "\n")
        else:
            tsv.write("\t".join(result) + "\n")
