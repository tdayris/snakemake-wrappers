#!/usr/bin/python3.8
# -*- coding: utf-8 -*-

"""
Annotate a VCF with several TSV from PharmkDB
"""


import datetime
import logging
import pandas
import numpy
import re

from typing import Any

logging.basicConfig(
    #filename=snakemake.log[0],
    filemode="w",
    level=logging.DEBUG
)


def get_headers(description: dict[str, Any]) -> str:
    """
    From a list of column name, and an optional list of description,
    build VCF headers.
    """
    version = 1.0
    name = "pharmkdb_annotate"
    url = f"github.com/tdayris/snakemake-wrappers/tree/Unofficial/bio/BiGR/{name}/wrapper.py"
    headers = [
        f"""##BiGRCommandLine=<ID={name},CommandLine="{url}",Version={version},Date={datetime.date.today()}>\n"""
    ]
    headers += [
        f"""##INFO<ID=f"PharmKDB_{key}",Number={value["nb"]},Type={value["type"]},Description="{value["desc"]}">\n"""
        for key, value in description.items()
    ]
    headers.append(
        f"""##FILTER<ID=ExistsInPharmkDBAnnMetadata,Description="Variant exists in PharmK DB Annotation Metadata">\n"""
    )
    return "".join(headers)


def dict_to_info(annot: dict[str, Any]) -> str:
    """Convert dict to VCF INFO annotation tag"""
    res = []
    for key, value in annot.items():
        if isinstance(value, bool):
            res.append(f"PharmKDB_{key}")
        elif (value == "") or (value is None):
            continue
        else:
            value = str(value).replace(";", ",").replace(", ", ",")
            res.append(f"PharmKDB_{key}={value}")

    return ";".join(res)


def annotate_metadata(line: str, tsv: pandas.DataFrame) -> str:
    """Annotate a variant from metadata tsv"""
    rs_id = line[:-1].split("\t")[2]

    if rs_id in [".", ""]:
        return line

    try:
        annot = dict_to_info(tsv.loc[rs_id].to_dict())
        chomp = line[:-1].split("\t")
        if chomp[7] == "":
            chomp[7] = annot
        else:
            chomp[7] += f";{annot}"

        if chomp[6] in [".", "", "PASS"]:
            chomp[6] = "ExistsInPharmkDBAnnMetadata"
        else:
            chomp[6] += ";ExistsInPharmkDBAnnMetadata"
        line = "\t".join(chomp) + "\n"
    except KeyError:
        logging.warning(f"No annotation for {rs_id} from {line}")

    return line


# Load CanceGeneCensus TSV file
logging.info("Loading annotation DB")
tsv  = pandas.read_csv(
    snakemake.input["clinical_ann_metadata"],
    sep="\t",
    header=0,
    index_col=None
)
logging.debug("Complete table, non modified")
logging.debug(tsv.head())

logging.debug("Location:")
logging.debug(tsv.Location.head())

logging.debug("Header:")
logging.debug(tsv.columns.tolist())


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
tsv.set_index("Location", inplace=True)
logging.debug(tsv.head())


description = {
    "Clinical_Annotation_ID": {
        "type": "Integer",
        "nb": "1",
        "desc": "Variant ID in PharmkDB"
    },
    "Location": {
        "type": "String",
        "nb": "1",
        "desc": "Variant rs ID in PharmkDB"
    },
    "Gene": {
        "type": "String",
        "nb": "1",
        "desc": "Gene name in PharmkDB"
    },
    "Level_of_Evidence": {
        "type": "String",
        "nb": "1",
        "desc": "Variant level of evidence in PharmkDB"
    },
    "Clinical_Annotation_Types": {
        "type": "String",
        "nb": "1",
        "desc": "Clinical annotation type in PharmkDB"
    },
    "Genotype_Phenotype_IDs": {
        "type": "String",
        "nb": "1",
        "desc": "Genotypes/Phenotypes described in PharmkDB"
    },
    "Annotation_Text": {
        "type": "String",
        "nb": "1",
        "desc": "Variant annotation and description in PharmkDB"
    },
    "Variant_Annotations_IDs": {
        "type": "String",
        "nb": "1",
        "desc": "Variant annotation ID in PharmkDB"
    },
    "Variant_Annotations": {
        "type": "String",
        "nb": "1",
        "desc": "Variant annotation in PharmkDB"
    },
    "PMIDs": {
        "type": "String",
        "nb": "1",
        "desc": "Pubmed IDs from PharmkDB"
    },
    "Evidence_Count": {
        "type": "String",
        "nb": "1",
        "desc": "Number of evidence as stated by PharmkDB"
    },
    "Related_Chemicals": {
        "type": "String",
        "nb": "1",
        "desc": "PharmkDB interaction between variant and chemicals"
    },
    "Related_Diseases": {
        "type": "String",
        "nb": "1",
        "desc": "PharmkDB interaction between variant and known diseases"
    },
    "Biogeographical_Groups": {
        "type": "String",
        "nb": "1",
        "desc": "Ethnicity on which studies were based, if any, as stated in PharmkDB"
    },
    "Chromosome": {
        "type": "String",
        "nb": "1",
        "desc": "PharmkDB Chromosome corresponding to variant"
    },
    "Latest_History_Date_YYYY_MM_DD": {
        "type": "String",
        "nb": "1",
        "desc": "PharmkDB last edition of variant description"
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
            new_headers = get_headers(description)
            line = new_headers + line
            logging.info("Header modified")

        else:
            line = annotate_metadata(line, tsv)

        out_vcf.write(line)
