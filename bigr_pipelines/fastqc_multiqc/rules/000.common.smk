# This file contains pure python functions only

import os
import re
import sys
from pathlib import Path
from typing import List
from snakemake.utils import min_version
min_version("7.0")


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

from snakemake.utils import min_version
from snakemake.shell import shell
min_version("6.0")

logging.basicConfig(
    filename="snakemake.fastqc_multiqc.log",
    filemode="w",
    level=logging.DEBUG
)

default_config = read_yaml(workflow_source_dir / "config" / "config.yaml")
config_path = get_config(default_config)
design = get_design(os.getcwd(), search_fastq_files)
try:
    design.columns = ["Sample_id", "Upstream_file"]
except ValueError:
    pass

fastq_links = link_fq(
    design.Sample_id,
    design.Upstream_file,
)

configfile: config_path
container: config.get("miniconda", "/mnt/beegfs/pipelines/unofficial-snakemake-wrappers/singularity/mambaforge_4.14.0-0.sif")


def get_archives(file: str, prefix: str = ".") -> List[str]:
    """
    Search file in prefix, recursively
    """
    print(f"Looking for {file} in {prefix}")
    prefix = Path(prefix)
    for content in prefix.iterdir():
        if content.is_dir() and content.name not in [".snakemake", "logs", "tmp", "data_output"]:
            yield from get_archives(file, str(content.absolute()))
        elif content.name == file:
            file_path = str(content.absolute())
            logging.info("%s was found", file_path)
            yield file_path

default_existing_path = config["miniconda"]
stats = next(
    get_archives(file="Stats.json.zip"),
    default_existing_path
)

interop = next(
    get_archives(file="InterOp.zip"),
    default_existing_path
)

runinfo = next(
    get_archives(file="RunInfo.xml.zip"),
    default_existing_path
)

runparams = next(
    get_archives(file="RunParameters.xml.zip"),
    default_existing_path
)

prefix = (
    "output"
    if any(zipfile != default_existing_path for zipfile in [interop, runinfo, runparams, stats])
    else "data_output"
)

def demux_input(stats, interop, runinfo, runparams) -> Dict[str, str]:
    """Return list of dmux-related expected input"""
    base = {
        "fqc_zip": expand(
            "fastqc/{sample}_fastqc.zip",
            sample=design["Sample_id"],
        ),
        "fqc_html": expand(
            "fastqc/{sample}.html",
            sample=design["Sample_id"],
        ),
        "txt": expand(
            "fastq_screen/{sample}.fastq_screen.txt",
            sample=design["Sample_id"],
        ),
        "png": expand(
            "fastq_screen/{sample}.fastq_screen.png",
            sample=design["Sample_id"],
        ),
    }
    if stats != default_existing_path:
        base["bcl_json"] = "tmp/Stats.json"
    else:
        logging.warning("Stats.json.zip was not found")

    if interop != default_existing_path:
        base["interop"] = "tmp/InterOp"
    else:
        logging.warning("InterOp.zip was not found")

    if runinfo != default_existing_path:
        base["runinfo"] = "tmp/RunInfo.xml"
    else:
        logging.warning("RunInfo.xml.zip was not found")
    
    if runparams != default_existing_path:
        base["runparams"] = "tmp/RunParameters.xml"
    else:
        logging.warning("RunParameters.xml.zip was not found")

    return base

