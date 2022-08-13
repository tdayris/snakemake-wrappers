#!/usr/bin/python3.8
# -*- coding: utf-8 -*-

"""
split ANN field in the INFO column in a VCF file
"""

import datetime
import gzip
import logging
import re

logging.basicConfig(
    filename=snakemake.log[0],
    filemode="w",
    level=logging.INFO
)


def open_function(file: str):
    """Return the correct opening function"""
    if file.endswith(".gz"):
        return gzip.open(file, "rb")
    return open(file, "r")


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


def sort_headers(lines) -> str:
    infos = []
    formats = []
    filters = []
    contigs = []
    fileformats = []
    others = []

    for line in lines:
        if line.startswith("##INFO"):
            infos.append(line)
        elif line.startswith("##FORMAT"):
            formats.append(line)
        elif line.startswith("##FILTER"):
            filters.append(line)
        elif line.startswith("##contig"):
            contigs.append(line)
        elif line.startswith("##fileformat"):
            fileformats.append(line)
        else:
            others.append(line)

    return "".join(fileformats + filters + formats + infos + others + contigs)

default_chr = list(map(str, range(23))) + ["MT", "X", "Y"]
if "default_chr" in snakemake.params.keys():
    default_chr = snakemake.params["default_chr"]

sample_rename_dict = snakemake.params.get("rename_sample", {})

version = 1.0
name = "fix_vcf"
url = f"github.com/tdayris/snakemake-wrappers/tree/Unofficial/bio/BiGR/{name}/wrapper.py"
headers = f"""##BiGRCommandLine=<ID={name},CommandLine="{url}",Version={version},Date={datetime.date.today()}>\n"""

header_list = []


logging.debug("Opening VCFs")
if str(snakemake.output["vcf"]).endswith("vcf.gz"):
    out_vcf = snakemake.output["vcf"][:-3]
else:
    out_vcf = snakemake.output["vcf"]

with (open_function(snakemake.input["vcf"]) as vcfin,
      open(out_vcf, "w", encoding="utf-8") as vcfout):
    for line in vcfin:
        if isinstance(line, bytes):
            line = line.decode("utf-8")

        if line.startswith("##"):
            # Header/formats/filters...
            #vcfout.write(line)
            if line[:-1] == '##FILTER=<ID=IsGermline,Number=.,Type=String,Description="Variant exists in Normal">':
                header_list.append('##FILTER=<ID=IsGermline,Description="Variant exists in Normal">\n')
            elif line[:-1] == '##FILTER=<ID=IsSomatic,Number=.,Type=String,Description="Variant does not exists in Normal, but exists in Tumor">':
                header_list.append('##FILTER=<ID=IsSomatic,Description="Variant does not exists in Normal, but exists in Tumor">\n')
            else:
                header_list.append(line)
            continue

        elif line.startswith("#"):
            # Column names
            vcfout.write(sort_headers(header_list))
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


if str(snakemake.output["vcf"]).endswith("vcf.gz"):
    logging.info(f"Compressing {out_vcf}")
    shell("pbgzip -c {out_vcf} > {snakemake.output['vcf']} 2> {log}")
    logging.info(f"Indexing {snakemake.output['call']}")
    shell("tabix -p vcf {snakemake.output['vcf']} >> {log} 2>&1")
    logging.info(f"Removing temporary file {out_vcf}")
    shell("rm --verbose {out_vcf} >> {log} 2>&1")
