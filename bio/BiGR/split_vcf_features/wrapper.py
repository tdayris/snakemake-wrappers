#!/usr/bin/python3.8
# -*- coding: utf-8 -*-

"""
split ANN field in the INFO column in a VCF file
"""

import datetime
import gzip
import logging

logging.basicConfig(
    filename=snakemake.log[0],
    filemode="w",
    level=logging.INFO
)

annotation_tag = snakemake.params.get("annotation_tag", "ANN=")


def open_function(file: str):
    """Return the correct opening function"""
    if file.endswith(".gz"):
        return gzip.open(file, "rb")
    return open(file, "r")

version = 1.0
name = "split_vcf_features"
url = f"github.com/tdayris/snakemake-wrappers/tree/Unofficial/bio/BiGR/{name}/wrapper.py"
header = f"""##BiGRCommandLine=<ID={name},CommandLine="{url}",Version="{version}",Date="{datetime.date.today()}">\n"""


if str(snakemake.output["call"]).endswith("vcf.gz"):
    out_vcf = snakemake.output["call"][:-3]
else:
    out_vcf = snakemake.output["call"]

with open(snakemake.input["call"], "r") as input_vcf, \
     open(out_vcf, "w") as splitted_vcf:

    # Copy headers without any modification
    for nb, line in enumerate(input_vcf):
        logging.debug(f"Line:\n{line}")
        if line.startswith("##"):
            splitted_vcf.write(line)
            continue

        if line.startswith("#"):
            splitted_vcf.write(header)
            splitted_vcf.write(line)
            continue

        # Split line in order to get info field
        chomp = line[:-1].split("\t")
        logging.debug(f"Line splitted:\n{chomp}")

        # Let user have an idea of the overall progress
        if nb % 10000 == 0:
            logging.info(
                f"{nb} lines processed, working on {chomp[0]}:{chomp[1]}"
            )

        # Get info field and extract both 'ann' and its position
        info_field = chomp[7].split(";")
        logging.debug(f"Info splitted:\n{info_field}")
        ann = None
        pos = 0
        for position, info in enumerate(info_field):
            if info.startswith(annotation_tag):
                ann = info[len(annotation_tag):]
                pos = position
                break
        else:
            # Handle case when no ANN is provided
            splitted_vcf.write(line)
            continue
        logging.debug(f"Ann = {ann}, position {pos}")

        # Break down annotation if and only if it contains multiple annotations
        if "," in ann:
            ann = ann.split(",")
            for feature in ann:
                logging.debug(f"working on {feature}")
                tmp = info_field.copy()
                tmp[pos] = f"{annotation_tag}{feature}"
                tmp = ";".join(tmp)
                chomp[7] = tmp
                splitted_vcf.write("\t".join(chomp) + "\n")
        else:
            # Handling case annotation contains only one feature
            splitted_vcf.write(line)


if str(snakemake.output["call"]).endswith("vcf.gz"):
    logging.info(f"Compressing {out_vcf}")
    shell("pbgzip -c {out_vcf} > {snakemake.output['call']} 2> {log}")
    logging.info(f"Indexing {snakemake.output['call']}")
    shell("tabix -p vcf {snakemake.output['call']} >> {log} 2>&1")
    logging.info(f"Removing temporary file {out_vcf}")
    shell("rm --verbose {out_vcf} >> {log} 2>&1")