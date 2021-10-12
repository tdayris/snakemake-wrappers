#!/usr/bin/python3.8
# -*- coding: utf-8 -*-

"""
Annotate a VCF with a TSV file from Cancer Gene Census
"""

import datetime
import logging
import pandas
import numpy
import re

from typing import Any

logging.basicConfig(
    filename=snakemake.log[0],
    filemode="w",
    level=logging.DEBUG
)


def get_headers(cols: list[str], description: dict[str, Any]) -> str:
    """
    From a list of column name, and an optional list of description,
    build VCF headers.
    """
    colnames = None
    version = 1.0
    name = "cancer_gene_census_annotate"
    url = f"github.com/tdayris/snakemake-wrappers/tree/Unofficial/bio/BiGR/{name}/wrapper.py"
    headers = [
        f"""##BiGRCommandLine=<ID={name},CommandLine="{url}",Version={version},Date={datetime.date.today()}>\n"""
    ]
    headers += [
        f"""##INFO<ID=f"CanceGeneCensus_{key}",Number={value["nb"]},Type={value["type"]},Description="{value["desc"]}">\n"""
        for key, value in description.items()
    ]
    headers.append(
        f"""##FILTER<ID=ExistsInCanceGeneCensus,Description="Transcript exists in Cancer Gene Census">\n"""
    )
    return "".join(headers)


def dict_to_info(annot: dict[str, Any]) -> str:
    """Convert dict to VCF INFO annotation tag"""
    res = []
    for key, value in annot.items():
        if isinstance(value, bool):
            res.append(f"CanceGeneCensus_{key}")
        elif (value == "") or (value is None):
            continue
        else:
            res.append(f"CanceGeneCensus_{key}={str(value)}")

    return ";".join(res)


def annotate(line: str, tsv: pandas.DataFrame) -> str:
    """
    Annotate a VCF formatted line with CanceGeneCensus
    """
    # Extract transcript if from VCF line (must be annotated with SnpEff!)
    try:
        gene_name = [
            i.split("=")[-1]
            for i in line[:-1].split("\t")[7].split(";")
            if i.startswith("ANN")
        ][0].split("|")[3].split(".")[0]
    except IndexError:
        logging.error(f"Could not find SnpEff annotation in: {line}")
    else:
        try:
            annot = dict_to_info(tsv.loc[gene_name].to_dict())
            chomp = line[:-1].split("\t")
            if chomp[7] == "":
                chomp[7] = annot
            else:
                chomp[7] += f";{annot}"

            if chomp[6] in [".", "", "PASS"]:
                chomp[6] = "ExistsInCanceGeneCensus"
            else:
                chomp[6] += ";ExistsInCanceGeneCensus"
            line = "\t".join(chomp) + "\n"
        except KeyError:
            logging.warning(f"No annotation for {gene_name} from {line}")
    return line


# Load CanceGeneCensus TSV file
logging.info("Loading annotation DB")
tsv  = pandas.read_csv(
    snakemake.input["cgc"],
    sep="\t",
    header=0,
    true_values=["yes", "Yes"]
)

# Some characters are not allowed in VCF. They are replaced here.
new_cols = []
for colname in tsv .columns.tolist():
    new_cols.append(colname.replace("-", "_")
                           .replace("/", "_")
                           .replace("\\", "_")
                           .replace(" ", "_")
                           .replace("(", "")
                           .replace(")", "")
                           .replace("#", "nb"))

tsv.columns = new_cols
tsv.set_index("Gene_Symbol", inplace=True)
logging.debug(tsv.head())

description = {
    "Gene_Symbol": {
        "type": "String",
        "nb": "1",
        "desc": "Gene symbol from Cance Gene Census"
    },
    "Name": {
        "type": "String",
        "nb": "1",
        "desc": "Gene name from Cancer Gene Census"
    },
    "Entrez_GeneId": {
        "type": "String",
        "nb": "1",
        "desc": "Entrez GeneId from Cancer Gene Census"
    },
    "Genome_Location": {
        "type": "String",
        "nb": "1",
        "desc": "Gene location in genome from Cancer Gene Census"
    },
    "Tier": {
        "type": "Integer",
        "nb": "1",
        "desc": "Cancer Gene Census Tier"
    },
    "Hallmark": {
        "type": "String",
        "nb": "1",
        "desc": "Cancer Gene Census indicates Gene belongs to Hallmark"
    },
    "Chr_Band": {
        "type": "String",
        "nb": "1",
        "desc": "Chromosome Band from Cancer Gene Census"
    },
    "Somatic": {
        "type": "String",
        "nb": "1",
        "desc": "Cancer Gene Census got this variant in somatic samples"
    },
    "Germline": {
        "type": "String",
        "nb": "1",
        "desc": "Cancer Gene Census got this variant in germline samples"
    },
    "Tumour_TypesSomatic": {
        "type": "String",
        "nb": "1",
        "desc": "Cancer Gene Census Tumor type"
    },
}


# Annotating input VCF
logging.debug("Opening VCFs")
with (open(snakemake.input["vcf"], "r") as in_vcf,
      open(snakemake.output["vcf"], 'w') as out_vcf):
    for line in in_vcf:
        if line.startswith("##"):
            pass

        elif line.startswith("#"):
            new_headers = get_headers(
                tsv .columns.tolist(),
                description
            )
            line = new_headers + line
            logging.info("Header modified")

        else:
            line = annotate(line, tsv)

        out_vcf.write(line)
