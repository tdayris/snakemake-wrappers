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
from yaml import dump, safe_load
from typing import Any, Dict, List


min_version("7.5")

import sys

# My own libraries
common = str(workflow_source_dir.absolute() / ".." / "common" / "python")
sys.path.append(common)

from dataframes import *
from file_manager import *
from files_linker import *
from reservation import *
from messages import message

with open("design.yaml", "r") as yamlstreamin:
    design = safe_load(yamlstreamin)

io_dict = dict(zip(design["output"], design["input"]))
irods_extra = design.get("irods_extra", config.get("irods_extra", "-vKrf"))
max_threads = design.get("threads", config.get("max_threads", 5))

wildcard_constraints:
    output_file=r"|".join(design["output"])