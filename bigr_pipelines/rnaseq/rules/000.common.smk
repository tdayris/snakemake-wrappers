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


def expected_targets(steps: Dict[str, Any]) -> Dict[str, Any]:
    """Return the list of expected output files"""
    results = {}
    if steps.get("qc", True) is True:
        results["qc"] = "data_output/multiqc/MultiQC.QC.html"

    if steps.get("quant", False) is True:
        results["quant"] = "data_output/multiqc/MultiQC.Salmon.html"

    if steps.get("dge", False) is True:
        results["dge"] = expand(
            "data_output/DEseq2/{comparison}/MultiQC.DEseq2.html",
            comparison=output_prefixes,
        )

    if steps.get("gsea", False) is True:
        results["clusterprofiler"] = expand(
            "data_output/GSEA/{comparison}/MultiQC.GSEA.html",
            comparison=output_prefixes,
        )

    if steps.get("immunedeconv", False) is True:
        results["immunedeconv"] = "data_output/MultiQC/ImmuneDeconv.html"

    if steps.get("fusions", False) is True:
        results["fusions"] = "data_output/multiqc/MultiQC.Star.Chimera.html"

    return results


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
fastq_links = link_fq(
    sample_names=design.Sample_id,
    r1_paths=design.Upstream_file,
    r2_paths=design.Downstream_file,
    prefix="data_input",
)

# A list that holds all comparisons made in DESeq2.
# This is done in order to avoid checkpoints
# This list contains: [(factor, test, ref), ...]
logging.info("Building DESeq2 globals...")
comparison_levels = list(
    yield_comps(
        complete_design=design,
        aggregate=config["deseq2"]["design"].get("columns_to_aggregate"),
        remove=config["deseq2"]["design"].get("columns_to_ignore"),
        contains=config["deseq2"]["design"].get("include_only"),
    )
)

factor_list=map(str, [i[0] for i in comparison_levels]),

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
if config.get("write_comparisons", True):
    with open("list_of_possible_comparisons.txt", "w") as outcomplist:
        outcomplist.write("\n".join(output_prefixes) + "\n")

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

# QC: PCAExplorer
# List of PCA axes to consider
axes_list = ["ax_1_ax_2", "ax_2_ax_3", "ax_3_ax_4"]
# Draw elipses on PCA ... or not, or both
elipsis_list = ["with_elipse", "without_elipse"]
# Up/down stream reads, this pipeline takes only pair-ended libraries
streams = ["1", "2"]
# Features types
features = ["genes", "transcripts"]


# Specific tools extensions
fastp_ext = ["html", "json"]
fastqscreen_ext = ["txt", "png"]


# Star Mappings
# Type of mapping performed by the pipepline
maptypes = ["variants", "chimera"]


# DEseq2 results content list
content_list = ["SortedOnLogFC", "SortedOnPadj", "Complete"]


# Immune deconv tool list
tool_list = ["cibersort", "cibersort_abs", "mcp_counter", "epic", "quantiseq", "xcell"]


# VDJ alignment: MIXcr
# List of possible segment export in mixcr
segment_export_list = ["vUsage", "jUsage", "isotypeUsage", "vjUsage"]


# GSEA: ClusterProfiler
# List of PPI databases for clusterProfiler
ppi_list = config["clusterprofiler"].get("ppi").keys()
# List of Gene oriented databases for clusterProfiler
gmt_list = config["clusterprofiler"].get("gmt").keys()
# List of gene set analysis methods
gse_method_list = ["enrich"]
# List of clusterprofiler plots
cprof_plots = ["barplot", "dotplot", "upsetplot"]
# List of possible key types
keytypes = ["ENSEMBL", "ENTREZID", "SYMBOL", "ENSEMBLPROT"]
db_key_dict = db_keytype(
    gmts=config["clusterprofiler"]["gmt"],
    ppis=config["clusterprofiler"]["ppi"],
)
database_keytypes_list = [f"{db}.{key}" for db, key in db_key_dict.items()]


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
    content=r"|".join(content_list),
    tool=r"|".join(tool_list),
    segment=r"|".join(segment_export_list),
    ppi=r"|".join(ppi_list),
    gmt=r"|".join(gmt_list),
    gse_method=r"|".join(gse_method_list),
    cprof_plot=r"|".join(cprof_plots),
    keytype=r"|".join(keytypes),
    ext=r"|".join(fastqscreen_ext + fastp_ext),
