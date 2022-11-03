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
from typing import Any, Dict, List

min_version("7.5")

import sys

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
from gmt import *
from messages import message


logging.basicConfig(
    filename="snakemake.variant_calling_ampliseq.log",
    filemode="w",
    level=logging.DEBUG
)

localrules: bigr_copy_fq
ruleorder: bwa_mem > bwa_fixmate_meta_bwa_mem

default_config = read_yaml(workflow_source_dir / "config.hg38.yaml")
configfile: get_config(default_config)
design = get_design(os.getcwd(), search_fastq_pairs)


wildcard_constraints:
    sample = r"|".join(design["Sample_id"]),
    stream = r"1|2|R1|R2"
    
    
fastq_links = link_fq(
    design.Sample_id,
    design.Upstream_file,
    design.Downstream_file
)