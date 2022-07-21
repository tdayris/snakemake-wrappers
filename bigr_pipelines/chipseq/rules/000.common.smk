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
worflow_source_dir = Path(snakemake.workflow.srcdir(".."))
common = str(worflow_source_dir / ".." / "common" / "python")
sys.path.append(common)

from dataframes import *
from file_manager import *
from files_linker import *
from graphics import *
from write_yaml import *
from messages import message

#####################
# Setup environment #
#####################

# Save output stream in a file
logging.basicConfig(filename="snakemake.rnaseq.log", filemode="w", level=logging.DEBUG)
logging.info("Additional utils loaded")


# Find and load configfile
default_config = read_yaml(worflow_source_dir / "config.hg38.yaml")


configfile: get_config(default_config=default_config)


logging.info("Config file loaded")


# Load design file and duplicate sample id as row name
design = get_design(dirpath=os.getcwd(), search_func=search_fastq_pairs)
design.set_index("Sample_id", inplace=True)
design["Sample_id"] = design.index.tolist()
logging.info("Design file loaded")

##################################
# Setup globals and fix wilcards #
##################################

# Links fastq paths provided by users and fastq paths used in this pipeline
# this is done in order to handle iRODS paths.
logging.info("Building globals...")
fastq_links = link_fq(
    sample_names=design.Sample_id,
    r1_paths=design.Upstream_file,
    r2_paths=design.Downstream_file,
    prefix="data_input",
)

sample_list = design.index.tolist()
peak_types = ["broadPeak", "narrowPeak"]
streams = ["1", "2"]


wildcard_constraints:
    sample=r"|".join(sample_list),
    stream=r"|".join(streams),
    peaktype=r"|".join(peak_types),


############################
### Resource reservation ###
############################

# Memory and time reservation
def get_resources_per_gb(
    wildcards, input, attempt, multiplier: int = 1, base: int = 0
) -> int:
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
        # ((input.size // 1_000_000_000) * attempt * multiplier) + base,
        # Case there is less than 1gb in input
        (multiplier * attempt) + base,
        1,
    )


logging.info("Preparing memory calls...")
# Explicit time reservations
get_45min_per_attempt = functools.partial(get_resources_per_gb, multiplier=45)
get_1h_per_attempt = functools.partial(get_resources_per_gb, multiplier=60)
get_2h_per_attempt = functools.partial(get_resources_per_gb, multiplier=60 * 2)
get_4h_per_attempt = functools.partial(get_resources_per_gb, multiplier=60 * 4)

# Explicit memory reservation
get_1gb_per_attempt = functools.partial(get_resources_per_gb, multiplier=1024)
get_2gb_per_attempt = functools.partial(get_resources_per_gb, multiplier=1024 * 2)
get_4gb_per_attempt = functools.partial(get_resources_per_gb, multiplier=1024 * 4)
get_8gb_per_attempt = functools.partial(get_resources_per_gb, multiplier=1024 * 8)
get_10gb_per_attempt = functools.partial(get_resources_per_gb, multiplier=1024 * 10)
get_75gb_and_5gb_per_attempt = functools.partial(
    get_resources_per_gb, multiplier=1024 * 5, base=1024 * 75
)
