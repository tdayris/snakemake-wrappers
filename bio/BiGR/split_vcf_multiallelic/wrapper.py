#!/usr/bin/python3.8
# -*- coding: utf-8 -*-

"""
split multiallelic sites into multiple lines
"""

import datetime
import logging

logging.basicConfig(
    filename=snakemake.log[0],
    filemode="w",
    level=logging.INFO
)

def split_formats(format_header,
                  allele_nb,
                  sample_names,
                  *format_content):
    """
    Split an unknown number of multiallelic-formats in corresponding number
    of lines

    In this case, we also handle the fact that there might be more than
    one sample and/or more than one caller in the VCF format
    """

    splitted = [
        [[] for _ in range(len(sample_names))]  # One for each sample/caller
        for _ in range(allele_nb)  # One for each allele
    ]
    format_iterator = list(zip(  # Must be a list, not an iterator
        format_header.split(":"),  # The header information
        *[i.split(":") for i in format_content]  # Per sample/caller information
    ))

    for sample_idx, sample_name in enumerate(sample_names):
        # Iterating in that order makes things easier for format speading
        for field_name, *field_values in format_iterator:
            current_field_values = field_values[sample_idx].split(",")
            if len(current_field_values) == allele_nb:
                # There is a format value for each allele
                for value_idx, value in enumerate(current_field_values):
                    splitted[value_idx][sample_idx] += [value]
            else:
                # There is a value for the multiallelic site
                for idx in range(allele_nb):
                    splitted[idx][sample_idx] += [field_values[sample_idx]]

    # We have to join the format fields before getting the values back!
    for kidx, k in enumerate(splitted):
        for lidx, l in enumerate(k):
            splitted[kidx][lidx] = ":".join(l)

    return splitted


def split_info(info: str, allele_nb: int) -> list[str]:
    """
    Split multiallelic info field
    """
    splitted = [[] for _ in range(allele_nb)]
    for info_field in info.split(";"):
        try:
            # The field has a key=val structure
            name, values = info_field.split("=")
        except ValueError:
            # The field has no name, it is a simple flag
            name = ""
            values = info_field

        # Remove empty fields
        # values = list(filter(lambda x: x != "", values.split(",")))
        values = [i for i in values.split(",") if i != ""]
        if len(values) == allele_nb:
            # The field has one value per allele
            for idx, value in enumerate(values):
                splitted[idx] += [f"{name}={value}" if name != "" else value]
        else:
            # The field has one value for all alleles
            for idx in range(allele_nb):
                splitted[idx] += [
                    f"{name}={';'.join(values)}"
                    if name != ""
                    else ";".join(values)
                ]

    return splitted


colnames = None
version = 1.0
name = "split_vcf_multiallelic"
url = f"github.com/tdayris/snakemake-wrappers/tree/Unofficial/bio/BiGR/{name}/wrapper.py"
header = f'##BiGRCommandLine=<ID={name},CommandLine="{url}",Version="{version}",Date="{datetime.date.today()}">\n'

with open(snakemake.input["vcf"]) as invcf, \
     open(snakemake.output["vcf"], 'w') as outvcf:
    for line in invcf:
        if line.startswith("##"):
            # Header content
            outvcf.write(line)
            logging.info("Annotation field processed")
            continue

        elif line.startswith("#"):
            # Column titles
            colnames = line[1:-1].split("\t")
            outvcf.write(header)
            outvcf.write(line)
            logging.info("Header processed")
            continue

        # Case line is a variant
        variant = {
           colname: content
           for colname, content in zip(colnames, line[:-1].split("\t"))
        }

        if "," in variant["ALT"]:
            # Multi-allelic position
            alleles = variant["ALT"].split(",")
            infos = split_info(variant["INFO"], len(alleles))
            sample_col_names = list(variant.keys())[9:]
            print(sample_col_names)
            formats = split_formats(variant["FORMAT"], len(alleles), colnames[9:], *[variant[i] for i in sample_col_names])
            if variant["FILTER"] == ".":
                filter = "BiGRSplittedMultiallelicSite"
            else:
                filter = variant["FILTER"] + ",BiGRSplittedMultiallelicSite"
            for alt, info, format in zip(alleles, infos, formats):
                outvcf.write("\t".join([
                    variant["CHROM"],
                    variant["POS"],
                    variant["ID"],
                    variant["REF"],
                    alt,
                    variant["QUAL"],
                    variant["FILTER"],
                    ";".join(info),
                    variant["FORMAT"],
                    *format
                ]) + "\n")

            continue

        # Mono allelic position
        outvcf.write(line)