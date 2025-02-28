#!/usr/bin/env python3
# coding: utf-8

import pandas

maf = pandas.read_csv(snakemake@input["maf"], sep="\t", header=0, index_col=None)

padding = snakemake.params.get("padding", 0)

chrom = snakemake.params.get("chrom", "Chromosome")
start = snakemake.params.get("start", "Start_Position")
end = snakemake.params.get("end", "End_Position")
bed_cols = [chrom, start, end]

bed = maf[bed_cols]
bed[start] -= padding
bed[end] += padding

if "bed3" in snakemake.output.keys():
    bed.to_csv(
        snakemake.output["bed3"],
        sep="\t",
        header=False,
        index=False
    )

bed_cols.append(snakemake.params.get("gene", "Hugo_Symbol"))
if "bed4" in snakemake.output.keys():
    bed = maf[bed_cols]
    bed[start] -= padding
    bed[end] += padding
    bed.to_csv(
        snakemake.output["bed4"],
        sep="\t",
        header=False,
        index=False
    )


bed_cols.append(snakemake.params.get("score", "Score"))
bed_cols.append(snakemake.params.get("strand", "Strand"))

if "bed6" in snakemake.output.keys():
    bed = maf[bed_cols]
    bed[start] -= padding
    bed[end] += padding
    bed.to_csv(
        snakemake.output["bed6"],
        sep="\t",
        header=False,
        index=False
    )

if "bed9" in snakemake.output.keys():
    bed = maf[bed_cols]
    bed[start] -= padding
    bed[end] += padding
    colors = {
        "MODIFIER": "O,153,0",
        "LOW": "0,0,255",
        "MODERATE": "255,128,0",
        "HIGH": "255,0,0",
    }
    impact = snakemake.params.get("impact", "IMPACT")
    bed["thickStart"] = bed[start] + padding
    bed["thickEnd"] = bed[end] - padding
    bed["itemRgb"] = [colors.get(i, "0,0,0") for i in maf[impact]]
    bed.to_csv(
        snakemake.output["bed9"],
        sep="\t",
        header=False,
        index=False
    )