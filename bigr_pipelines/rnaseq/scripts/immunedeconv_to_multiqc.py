#!/usr/bin/env python3
# coding: utf-8

"""
This script takes TSV file from ImmundeDeconv
and returns a configuration file for MultiQC
"""


import logging
import pandas
import yaml

from typing import Any, Dict, Optional


def read_immunedeconv(path: str) -> pandas.DataFrame:
    """Load immunedeconv in memory"""
    logging.info(f"Loading dataframe at: {path}")
    # Load data frame
    tmp = pandas.read_csv(
        path,
        sep="\t",
        header=0,
        index_col=0,
    )

    # Drop lines full of zeros
    tmp = tmp.loc[~(tmp==0).all(axis=1)]
    logging.debug(tmp.head())
    return tmp


def load_config(path: str) -> Dict[str, Any]:
    """Load an existing MultiQC configuration file"""
    logging.info(f"Loading config at: {path}")
    with open(path, "r") as yaml_stream:
        yaml_data = yaml.safe_load(stream=yaml_stream)
        logging.debug(yaml_data)
        return yaml_data



def get_permissions(tool_name: str) -> str:
    """
    From the name of a tool handled by ImmuneDeconv, 
    return comparison permissions
    """
    logging.info(f"Building permission text for {tool_name}:")
    text = ""
    between_samples = [
        "mcp_counter", 
        "xcell", 
        "timer", 
        "consensus_tme", 
        "estimate", 
        "absis", 
        "mmcp_counter", 
        "base"
    ]

    between_cell_types = ["cibersort", "dcq"]

    both = ["epic", "quantiseq", "cibersort_abs", "seqimmucc"]

    if tool_name.lower() in between_samples:
        text = (
            f"{tool_name} method allows between-sample comparisons only, "
            "not between-cell-types."
        )

    if tool_name.lower() in between_cell_types:
        text = (
            f"{tool_name} method allows between-cell-type comparisons, "
            "not between-sample."
        )

    if tool_name.lower() in both:
        text = (
            f"{tool_name} method allows both between-sample "
            "and between-cell-type comparisons."
        )

    logging.debug(text)
    return text


def get_ylab(tool_name: str) -> str:
    """Return correct y lab term"""
    return "fraction"


def get_description(tool_name: str, 
                    tumor: Optional[str] = None) -> str:
    """
    Return a description for a given tool, 
    on an optional given tumor type
    """
    logging.info(f"Description for {tool_name} ({tumor}):")
    text = f"This deconvolution has been made by {tool_name}. "
    if tumor:
        text += (
            "It was performed against a known "
            f"dataset composed of TCGA's {tumor}. "
        )
    
    text += get_permissions(tool_name)
    logging.debug(text)
    return text
    

def build_config(tool_name: str, 
                 deconv_path: str, 
                 tumor: Optional[str] = None) -> Dict[str, Any]:
    """
    From a tool, path to its results and tumor name, 
    return a custom config section
    """
    logging.info("Building config")
    data = read_immunedeconv(path=deconv_path).to_dict()
    logging.debug(data)
    config = {
        "id": tool_name,
        "description": get_description(tool_name=tool_name, tumor=tumor),
        "plot_type": "bargraph",
        "pconfig": {
            "id": f"{tool_name.lower()}_barplot",
            "title": f"Call-type deconvolution with {tool_name}",
            "ylab": get_ylab(tool_name=tool_name),
        },
        "data": {data}
    }
    logging.debug(config)
    return config


# Main program
logging.basicConfig(
    filename=snakemake.log[0],
    filemode="w",
    level=logging.DEBUG
)

tumor_abbreviation = {
    "kich": "Kidney Chromophobe",
    "blca": "Bladder Urothelial Carcinoma",
    "brca": "Breast invasive carcinoma",
    "cesc": "Cervical squamous cell carcinoma and endocervical adenocarcinoma",
    "gbm": "Glioblastoma multiforme",
    "hnsc": "Head and Neck squamous cell carcinoma",
    "kirp": "Kidney renal papillary cell carcinoma",
    "lgg": "Brain Lower Grade Glioma",
    "lihc": "Liver hepatocellular carcinoma",
    "luad": "Lung adenocarcinoma",
    "lusc": "Lung squamous cell carcinoma",
    "prad": "Prostate adenocarcinoma",
    "sarc": "Sarcoma",
    "pcpg": "Pheochromocytoma and Paraganglioma",
    "paad": "Pancreatic adenocarcinoma",
    "tgct": "Testicular Germ Cell Tumors",
    "ucec": "Uterine Corpus Endometrial Carcinoma",
    "ov": "Ovarian serous cystadenocarcinoma",
    "skcm": "Skin Cutaneous Melanoma",
    "dlbc": "Lymphoid Neoplasm Diffuse Large B-cell Lymphoma",
    "kirc": "Kidney renal clear cell carcinoma",
    "acc": "Adrenocortical carcinoma",
    "meso": "Mesothelioma",
    "thca": "Thyroid carcinoma",
    "uvm": "Uveal Melanoma",
    "ucs": "Uterine Carcinosarcoma",
    "thym": "Thymoma",
    "esca": "Esophageal carcinoma",
    "stad": "Stomach adenocarcinoma",
    "read": "Rectum adenocarcinoma",
    "coad": "Colon adenocarcinoma",
    "chol": "Cholangiocarcinoma",
}

prefix = snakemake.params.get("prefix", "data_output/")
suffix = snakemake.params.get("suffix", ".tsv")
multiqc_config = {
    "title": "Immune cell type deconvolution", 
    "custom_data": {},
    "show_analysis_paths": False,
    "show_analysis_time": False,
    "report_header_info": [
        {"Contact E-mail": "bigr@gustaveroussy.fr"},
        {"Application Type": "RNA-seq"},
        {"Project Type": "Immune cell type deconvolution"},
    ],
}

for deconv_result in snakemake.input:
    logging.debug(deconv_result)

    if deconv_result.startswith(prefix):
        names = deconv_result[len(prefix):-len(suffix)]

    tool_name, tumor = names.split("/")
    logging.info(f"Tool name: {tool_name}")
    logging.info(f"Tumor abbrv: {tumor}")

    tumor = tumor_abbreviation.get(tumor)
    logging.info(f"Tumor name: {tumor}")

    multiqc_config["custom_data"][f"{tool_name.lower()}_results"] = build_config(
        tool_name=tool_name, 
        deconv_path=deconv_result, 
        tumor=tumor
    )
    logging.info("Config built:")
    logging.debug(multiqc_config)

with open(snakemake.output["yaml"], "w") as yaml_outstream:
    yaml.safe_dump(
        data=multiqc_config, 
        stream=yaml_outstream, 
        default_flow_style=False
    )

logging.info("Process over")