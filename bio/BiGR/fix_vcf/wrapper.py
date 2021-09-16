#!/usr/bin/python3.8
# -*- coding: utf-8 -*-

"""
split ANN field in the INFO column in a VCF file
"""

import datetime
import logging
import re

logging.basicConfig(
    filename=snakemake.log[0],
    filemode="w",
    level=logging.INFO
)

def create_header(section: str,
                  key: str,
                  tp: str,
                  nb: str,
                  desc: str) -> str:
    return f'##{section}=<ID={key},Number={nb},Type={tp},Description="{desc}">'


def add_info(key: str, value: str, info: str) -> str:
    if info != ".":
        info += ";"
    else:
        info = ""

    if key in info:
        logging.warning(f"Adding existing key: {key} to {info}")

    return f"{info}{key}={value}"


def add_format_header(key: str, format: str) -> str:
    return ":".join(map(str, [format, key]))

def add_sample(value: str, format: str) -> str:
    return ":".join(map(str, [format, value]))

def add_format(format_header: str,
               key: str,
               value: str,
               *formats: list[str]) -> str:
    result = [add_format_header(key, format_header)]
    for format in formats:
        result.append(add_sample(value, format))
    return "\t".join(result)

default_chr = list(map(str, range(23))) + ["MT", "X", "Y"]
if "default_chr" in snakemake.params.keys():
    default_chr = snakemake.params["default_chr"]

sample_rename_dict = snakemake.params.get("rename_sample", {})

version = 1.0
name = "fix_vcf"
url = f"github.com/tdayris/snakemake-wrappers/tree/Unofficial/bio/BiGR/{name}/wrapper.py"
headers = f"""##BiGRCommandLine=<ID={name},CommandLine="{url}",Version={version},Date={datetime.date.today()}>\n"""


with (open(snakemake.input["vcf"], "r") as vcfin,
      open(snakemake.output["vcf"], "w") as vcfout):
    for line in vcfin:
        if line.startswith("##"):
            # Header/formats/filters...
            vcfout.write(line)
            continue

        elif line.startswith("#"):
            # Column names
            vcfout.write(headers)
            print(line)
            chomp = line[:-1].split("\t")
            print(chomp)
            for idx, col in enumerate(chomp):
                chomp[idx] = sample_rename_dict.get(col, col)
            line = "\t".join(chomp) + "\n"
            vcfout.write(line)
            continue

        chrom, pos, idx, ref, alt, qual, fil, info, format, *samples = line[:-1].split("\t")

        if (snakemake.params.get("remove_non_conventional_chromosomes", True) is True) and (chrom not in default_chr):
            continue

        info = info.replace(";;", ";").strip(";")
        line = "\t".join([chrom, pos, idx, ref, alt, qual, fil, info, format, *samples]) + "\n"
        vcfout.write(line)
