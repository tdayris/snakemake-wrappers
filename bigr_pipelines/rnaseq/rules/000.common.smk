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

# My own libraries
worflow_source_dir = Path(snakemake.workflow.srcdir("."))
common = str(worflow_source_dir / ".." / "common" / "python")
sys.path.append(common)

from file_manager import *
from files_linker import *
from write_yaml import *
from messages import message


#####################
# Setup environment #
#####################

# Save output stream in a file
logging.basicConfig(
    filename="snakemake.rnaseq.log", filemode="w", level=logging.DEBUG
)


# Find and load configfile
default_config = read_yaml(worflow_source_dir / "config.hg38.yaml")
configfile: get_config(default_config=default_config)

# Load design file
design = get_design(dirpath=os.getcwd(), search_func=search_fastq_pairs)

##################################
# Setup globals and fix wilcards #
##################################

# Links fastq paths provided by users and fastq paths used in this pipeline
# this is done in order to handle iRODS paths.
fastq_links = link_fq(design.Sample_id, design.Upstream_file, design.Downstream_file)

# Globals used in wildcards
sample_list = design.Sample_id.to_list()
streams = ["1", "2"]


wildcard_constraints:
    sample=r"|".join(sample_list),
    stream=r"|".join(streams),


############################
### Resource reservation ###
############################

# Memory and time reservation
def get_resources_per_gb(wildcards, input, attempt, multiplier) -> int:
    """
    Return the amount of resources needed per GB of input.

    Parameters:
    * wildcards: Snakemake signature requires this parameter.
    * input: The input file
    * attempt: The # of times the calling rule has been restarted
    * multiplier: An arbitrary multiplier

    Return:
    (int) The amount of resources needed (mb, minutes, etc)
    """
    return max(
        # Case there is 1gb or more in input
        (input.size // 1_000_000_000) * attempt * multiplier,
        # Case there is less than 1gb in input
        multiplier * attempt
    )

# Explicit time reservations
get_75min_per_attempt = functools.partial(get_resources_per_gb, multiplier=15)
get_45min_per_attempt = functools.partial(get_resources_per_gb, multiplier=45)
get_75min_per_attempt = functools.partial(get_resources_per_gb, multiplier=75)

# Explicit memory reservation
get_1gb_per_attempt = functools.partial(get_resources_per_gb, multiplier=1024)
get_4gb_per_attempt = functools.partial(get_resources_per_gb, multiplier=1024 * 4)
get_10gb_per_attempt = functools.partial(get_resources_per_gb, multiplier=1024 * 10)