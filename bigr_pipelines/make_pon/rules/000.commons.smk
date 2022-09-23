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
from reservation import *
from messages import message

#####################
# Setup environment #
#####################

# Save output stream in a file
logging.basicConfig(filename="snakemake.rnaseq.log", filemode="w", level=logging.DEBUG)
logging.info("Additional utils loaded")


# Find and load configfile
default_config = read_yaml(workflow_source_dir / "config" / "config.yaml")


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


genome_id = config.get("used_genome", "hg38")
fasta_file = config[genome_id]["fasta"]


# Given a fasta file: /path/to/sequence.fasta
# The fasta index should be: /path/to/sequence.fasta.fai
fai_file = config[genome_id]["fasta"] + ".fai"

# Given a fasta file: /path/to/sequence.fasta
# The sequence dictionary should be: /path/to/sequence.dict
dict_file = config[genome_id].get(
    "dict", ".".join(os.path.splitext(fasta_file)[:-1]) + ".dict"
)

# Given a vcf file: /path/to/known.vcf.gz
# The tabbix should be: /path/to/known.vcf.gz.tbi
tbi_file = config[genome_id].get("tbi", config[genome_id]["vcf"] + ".tbi")

# Given a fasta file: /path/to/sequence.fasta
# The bwa index may be: /path/to/sequence.{0123,amb,ann,bwt.2bit.64,pac}
bwa_sequence_index = config[genome_id].get(
    "bwa_index",
    multiext(
        ".".join(os.path.splitext(fasta_file)[:-1]),
        ".0123",
        ".amb",
        ".ann",
        ".bwt.2bit.64",
        ".pac",
    ),
)

# Loading list of samples used in the PoN
design = pandas.read_csv(
    config["design"], sep="\t", header=0, index_col=None
).set_index("Sample_id")
