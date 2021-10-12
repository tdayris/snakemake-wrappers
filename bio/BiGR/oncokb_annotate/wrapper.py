#!/usr/bin/python3.8
# -*- coding: utf-8 -*-

"""
Annotate a VCF with a CSV file from OncoKB
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
    name = "oncokb_annotate"
    url = f"github.com/tdayris/snakemake-wrappers/tree/Unofficial/bio/BiGR/{name}/wrapper.py"
    headers = [
        f"""##BiGRCommandLine=<ID={name},CommandLine="{url}",Version={version},Date={datetime.date.today()}>\n"""
    ]
    headers += [
        f"""##INFO<ID=f"OncoKB_{key}",Number={value["nb"]},Type={value["type"]},Description="{value["desc"]}">\n"""
        for key, value in description.items()
    ]
    headers.append(
        f"""##FILTER<ID=ExistsInOncoKB,Description="Transcript exists in OncoKB">\n"""
    )
    return "".join(headers)


def dict_to_info(annot: dict[str, Any]) -> str:
    """Convert dict to VCF INFO annotation tag"""
    res = []
    for key, value in annot.items():
        if isinstance(value, bool):
            res.append(f"OncoKB_{key}")
        elif (value == "") or (value is None):
            continue
        else:
            res.append(f"OncoKB_{key}={str(value)}")

    return ";".join(res)


def annotate(line: str, csv: pandas.DataFrame) -> str:
    """
    Annotate a VCF formatted line with OnkoKB
    """
    # Extract transcript if from VCF line (must be annotated with SnpEff!)
    try:
        ensembl_transcript = [
            i.split("=")[-1]
            for i in line[:-1].split("\t")[7].split(";")
            if i.startswith("ANN")
        ][0].split("|")[6].split(".")[0]
    except IndexError:
        logging.error(f"Could not find SnpEff annotation in: {line}")
    else:
        try:
            annot = dict_to_info(csv.loc[ensembl_transcript].to_dict())
            chomp = line[:-1].split("\t")
            if chomp[7] == "":
                chomp[7] = annot
            else:
                chomp[7] += f";{annot}"

            if chomp[6] in [".", "", "PASS"]:
                chomp[6] = "ExistsInOncoKB"
            else:
                chomp[6] += ";ExistsInOncoKB"
            line = "\t".join(chomp) + "\n"
        except KeyError:
            logging.warning(f"No annotation for {ensembl_transcript} from {line}")
            pass

    return line


# Load OncoKB CSV file
logging.info("Loading annotation DB")
csv  = pandas.read_csv(
    snakemake.input["oncokb"],
    sep=",",
    header=0,
    true_values=["Yes"],
    false_values=["No"]
)

# Some characters are not allowed in VCF. They are replaced here.
new_cols = []
for colname in csv .columns.tolist():
    new_cols.append(colname.replace("-", "_")
                           .replace("/", "_")
                           .replace("\\", "_")
                           .replace(" ", "_")
                           .replace("(", "")
                           .replace(")", "")
                           .replace("#", "nb"))

csv.columns = new_cols
csv.set_index("GRCh38_Isoform", inplace=True)

description = {
    "Hugo_Symbol": {
        "type": "String",
        "nb": "1",
        "desc": "Onco KB Hugo Symbol"
    },
    "Entrez_Gene_ID": {
        "type": "String",
        "nb": "1",
        "desc": "Onco KB Entrez Gene ID"
    },
    "GRCh37_Isoform": {
        "type": "String",
        "nb": "1",
        "desc": "Onco KB Ensembl GRCh37 transcript ID"
    },
    "GRCh37_RefSeq": {
        "type": "String",
        "nb": "1",
        "desc": "Onco KB Refseq ID in GRCh37"
    },
    "GRCh38_RefSeq": {
        "type": "String",
        "nb": "1",
        "desc": "Onco KB Refseq ID in GRCh38"
    },
    "GRCh37_Isoform": {
        "type": "String",
        "nb": "1",
        "desc": "Onco KB Ensembl GRCh38 transcript ID"
    },
    "#_of_occurrence_within_resources_Column_D_J": {
        "type": "Integer",
        "nb": "1",
        "desc": "Onco KB number of occurrence within resources Column DJ"
    },
    "OncoKB_Annotated": {
        "type": "Flag ",
        "nb": "0",
        "desc": "The gene is annotated by OncoKB"
    },
    "Is_Oncogene": {
        "type": "Flag ",
        "nb": "0",
        "desc": "The gene is annotated as an oncogene by OncoKB"
    },
    "Is_Tumor_Suppressor_Gene": {
        "type": "Flag ",
        "nb": "0",
        "desc": "The gene is annotated as tumor suppressor by OncoKB"
    },
    "MSK_IMPACT": {
        "type": "Flag ",
        "nb": "0",
        "desc": "Gene has OncoKB integrated mutation profiling of actionable cancer targets"
    },
    "MSK_HEME": {
        "type": "Flag ",
        "nb": "0",
        "desc": "Gene has OncoKB integrated mutation profiling of actionable cancer targets within blood"
    },
    "FOUNDATION_ONE": {
        "type": "Flag ",
        "nb": "0",
        "desc": "The gene is annotated by OncoKB as belonging to FoundationOne CDx"
    },
    "FOUNDATION_ONE_HEME": {
        "type": "Flag ",
        "nb": "0",
        "desc": "The gene is annotated by OncoKB as belonging to FoundationOne CDx, blood samples"
    },
    "Vogelstein": {
        "type": "Flag ",
        "nb": "0",
        "desc": "The gene is annotated by OncoKB as belonging Vogelstein 2013 publication"
    },
    "SANGER_CGC05_30_2017": {
        "type": "Flag ",
        "nb": "0",
        "desc": "The gene is annotated by OncoKB as belonging to Cancer Gene Sensus"
    }
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
                csv .columns.tolist(),
                description
            )
            line = new_headers + line
            logging.info("Header modified")

        else:
            line = annotate(line, csv)

        out_vcf.write(line)
