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
workflow_source_dir = Path(snakemake.workflow.srcdir(".."))
common = str(workflow_source_dir / ".." / "common" / "python")
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
default_config = read_yaml(workflow_source_dir / "config.hg38.yaml")


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
fastq_links = link_fq(design.Sample_id, design.Upstream_file, design.Downstream_file)

# A list that holds all comparisons made in DESeq2.
# This is done in order to avoid checkpoints
# This list contains: [(factor, test, ref), ...]
logging.info("Building DESeq2 globals...")
comparison_levels = list(
    yield_comps(
        complete_design=design,
        aggregate=config["deseq2"]["design"].get("aggregate_col"),
        remove=config["deseq2"]["design"].get("remove_col"),
        contains=config["deseq2"]["design"].get("include_only"),
    )
)

# An iterator that holds all samples involved in the comparisons
# listed above
samples_iterator = yield_samples(
    complete_design=design.copy(),
    aggregate=config["deseq2"]["design"].get("aggregate_col"),
    remove=config["deseq2"]["design"].get("remove_col"),
)

# A list containing all expected PCA at the end of DESeq2.
# This is done in order to avoid checkpoints
expected_pcas = [
    f"figures/DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}/pca/pca_{factor}_{axes}_{elipse}.png"
    for (factor, test, ref) in comparison_levels
    for axes in ["ax_1_ax_2", "ax_2_ax_3"]  # , "ax_3_ax_4"]
    for elipse in ["with_elipse", "without_elipse"]
]

# A dict containing sample levels for a given factor.
condition_dict = {
    f"DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}": relation_condition_sample(
        design.copy(), factor, test, ref
    )
    for factor, test, ref in comparison_levels
}

# A list containing all DESeq2 comparison names
# Stored as a list for futrther re-use
output_prefixes = [
    f"DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}"
    for factor, test, ref in comparison_levels
]

# A dict containing comparison names and factors/levels
contrasts = dict(zip(output_prefixes, comparison_levels))

# A dict containing the list of samples used in DESeq2,
# for a given comparison
samples_per_prefixes = {
    prefix: list(condition_dict[prefix].keys()) for prefix in output_prefixes
}

# Boolean: is there any batch effect to consider in DEseq2 ?
batch_effect = any(level[0] == "BatchEffect" for level in comparison_levels)


# Globals used in wildcards
# List of samples
sample_list = design.Sample_id.to_list()
# List of PCA axes to consider
axes_list = ["ax_1_ax_2", "ax_2_ax_3", "ax_3_ax_4"]
# Draw elipses on PCA ... or not, or both
elipsis_list = ["with_elipse", "without_elipse"]
# Up/down stream reads, this pipeline takes only pair-ended libraries
streams = ["1", "2"]
# Features types
features = ["gene", "transcript"]
# Type of mapping performed by the pipepline
maptypes = ["variants", "chimera"]

logging.info("Constraining wildcards...")


wildcard_constraints:
    sample=r"|".join(sample_list),
    stream=r"|".join(streams),
    comparison=r"|".join(output_prefixes),
    factor=r"|".join(map(str, [i[0] for i in comparison_levels])),
    test=r"|".join(map(str, [i[1] for i in comparison_levels])),
    ref=r"|".join(map(str, [i[2] for i in comparison_levels])),
    axes=r"|".join(axes_list),
    elipse=r"|".join(elipsis_list),
    feature=r"|".join(features),
    maptype=r"|".join(maptypes),


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
        ((input.size // 1_000_000_000) * attempt * multiplier) + base,
        # Case there is less than 1gb in input
        (multiplier * attempt) + base,
    )


logging.info("Preparing memory calls...")
# Explicit time reservations
get_15min_per_attempt = functools.partial(get_resources_per_gb, multiplier=15)
get_45min_per_attempt = functools.partial(get_resources_per_gb, multiplier=45)
get_75min_per_attempt = functools.partial(get_resources_per_gb, multiplier=75)
get_20min_per_attempt = functools.partial(get_resources_per_gb, multiplier=20)
get_1h_per_attempt = functools.partial(get_resources_per_gb, multiplier=60)
get_2h_per_attempt = functools.partial(get_resources_per_gb, multiplier=60 * 2)
get_5h_per_attempt = functools.partial(get_resources_per_gb, multiplier=60 * 5)
get_90min_per_attempt = functools.partial(get_resources_per_gb, multiplier=90)

# Explicit memory reservation
get_1gb_per_attempt = functools.partial(get_resources_per_gb, multiplier=1024)
get_2gb_per_attempt = functools.partial(get_resources_per_gb, multiplier=1024 * 2)
get_4gb_per_attempt = functools.partial(get_resources_per_gb, multiplier=1024 * 4)
get_6gb_per_attempt = functools.partial(get_resources_per_gb, multiplier=1024 * 6)
get_8gb_per_attempt = functools.partial(get_resources_per_gb, multiplier=1024 * 8)
get_10gb_per_attempt = functools.partial(get_resources_per_gb, multiplier=1024 * 10)
get_20gb_per_attempt = functools.partial(get_resources_per_gb, multiplier=1024 * 20)
get_75gb_and_5gb_per_attempt = functools.partial(
    get_resources_per_gb, multiplier=1024 * 5, base=1024 * 75
)
