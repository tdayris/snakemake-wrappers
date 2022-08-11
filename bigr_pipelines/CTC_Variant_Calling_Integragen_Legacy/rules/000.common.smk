"""
This snakefile contains python functions and globals
Skip if you're looking for rules
"""

#####################################
# Check Snakemake version           #
# Import search-and-buil functions  #
# for config and design             #
#####################################

# Official libraries
import os
import functools

from snakemake.utils import min_version
from pathlib import Path
from yaml import dump

min_version("7.5")

import sys
import pandas

# My own libraries
workflow_source_dir = Path(snakemake.workflow.srcdir(".."))
common = str(workflow_source_dir / ".." / "common" / "python")
sys.path.append(common)

from dataframes import *
from file_manager import *
from files_linker import *
from graphics import *
from write_yaml import *
from reservation import *
from messages import message

#####################
# Setup environment #
#####################

# Save output stream in a file
logging.basicConfig(filename="snakemake.rnaseq.log", filemode="w", level=logging.DEBUG)
logging.info("Additional utils loaded")


# Find and load configfile
default_config = read_yaml(workflow_source_dir / "config.hg38.yaml")


configfile: get_config(default_config=default_config)


logging.info("Config file loaded")


# Load design file and duplicate sample id as row name
design = get_design(dirpath=os.getcwd(), search_func=search_mapping)
design.set_index("Sample_id", inplace=True)
design["Sample_id"] = design.index.tolist()
logging.info("Design file loaded")

##################################
# Setup globals and fix wilcards #
##################################

# Links fastq paths provided by users and fastq paths used in this pipeline
# this is done in order to handle iRODS paths.
logging.info("Building globals...")


def parse_design(
    design: pandas.DataFrame, prefix: str = "data_input", suffix: str = "bam"
) -> dict[str, str]:

    link_bams = {}
    sample_list = []
    link_sample_baseline = {}

    row = next(design.iterrows(), None)[1]

    while row is not None:
        if row["Status"].lower() == "baseline":
            result[f"{prefix}/{sample}.baseline.{suffix}"] = row["bam"]

        elif row["Status"].lower() == "wbc":
            manip = row["Manip"]
            kit = row["Version"]
            sample_id = f"{sample}_V{kit}_M{manip}"

            result[f"{prefix}/{sample_id}.wbc.{suffix}"] = row["bam"]

        elif row["Status"].lower() == "ctc":
            manip = row["Manip"]
            kit = row["Version"]
            replicate = row["Replicate"]
            raw_sample_id = f"{sample}_V{kit}_M{manip}"
            sample_id = f"{raw_sample_id}_{replicate}"

            sample_list.append(sample_id)
            result[f"{prefix}/{sample_id}.ctc.{suffix}"] = row["bam"]
            link_sample_baseline[sample_id] = {
                "ctc": f"{prefix}/{sample_id}.ctc.{suffix}",
                "wbc": f"{prefix}/{raw_sample_id}.wbc.{suffix}",
                "baseline": f"{prefix}/{sample}.baseline.{suffix}",
            }

    return link_bams, sample_list, link_sample_baseline


def get_baseline(wildcards):
    return link_sample_baseline[wildcards.sample]["baseline"]


def get_wbc(wildcards):
    return link_sample_baseline[wildcards.sample]["wbc"]


def get_ctc(wildcards):
    return link_sample_baseline[wildcards.sample]["ctc"]


def get_trio(wildcards):
    return {
        "tumor": link_sample_baseline[wildcards.sample]["ctc"],
        "wbc": link_sample_baseline[wildcards.sample]["wbc"],
        "normal": link_sample_baseline[wildcards.sample]["baseline"],
        "fasta": config["ref"]["fasta"],
    }


link_bams, sample_list, link_sample_baseline = parse_design(design.copy())
replicate_list = list(set(design["Replicate"]))
version_list = list(set(design["Version"]))
manip_list = list(set(design["Manip"]))
status_list = list(set(design["Status"]))


wildcard_constraints:
    sample=r"|".join(samples_list),
    status=r"|".join(status_list),