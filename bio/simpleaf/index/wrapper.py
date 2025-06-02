# coding: utf-8

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2025, Thiault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

import os

from snakemake.shell import shell
from distutils.spawn import find_executable
from pathlib import Path
from snakemake.shell import shell
from tempfile import TemporaryDirectory
from warnings import warn

log = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)
extra = snakemake.params.get("extra", "")

# ALEVIN_FRY_HOME environment variable is required for
# simpleaf to run.
alevin_fry_home = os.environ.get("ALEVIN_FRY_HOME")
if alevin_fry_home:
    # Case it was provided
    alevin_fry_home = Path(alevin_fry_home)
else:
    # Case the environment variable is missing in shell env
    alevin_fry_home = snakemake.params.get(
        "alevin_fry_home", Path(snakemake.output[0]).parent
    )
    os.environ["ALEVIN_FRY_HOME"] = str(alevin_fry_home.absolute())

    # Tell used where the ALEVIN_FRY_HOME
    warn(
        "ALEVIN_FRY_HOME variable not found in environment. "
        f"It has been set to {os.environ['ALEVIN_FRY_HOME']}"
    )


# Simpleaf configuration file MUST be available in
# ALEVIN_FRY_HOME directory. Its name must be: simpleaf_info.json
info_path = snakemake.input.get("info")
if info_path:
    # Make info available in ALEVIN_FRY_HOME directory
    info_path = Path(info_path)
    if info_path.parent != alevin_fry_home:
        # Then info_path is not into alevin_fry_home and must be linked
        os.symlink(
            src=str(info_path.absolute()),
            dest=str(alevin_fry_home.absolute()),
        )
else:
    # No existing simpleaf_info.json, it must be created.
    # It will create the output file into $ALEVIN_FRY_HOME, which has been
    # defined above.
    shell("simpleaf set-paths {log}")

gtf = snakemake.input.get("gtf")
if gtf and str(gtf).lower().endswith((".gff", ".gff3")):
    extra += " --gff3-format "

index_input = ""
# User is expected to know mutually exclusive input files
# or simpleaf itself will raise a proper error.
expected_input = {
    "gtf": "--gtf",
    "feature": "--feature-csv",
    "probe": "--probe-csv",
    "refseq": "--ref-seq",
    "fasta": "--fasta",
    "spliced": "--spliced",
    "unspliced": "--unspliced",
    "decoy": "--decoy-paths",
    "info": None,
    "chemistry": None,
}
for snakemake_key, cmd_arg in expected_input.items():
    tmp = snakemake.input.get(snakemake_key)
    if tmp and (cmd_arg is not None):
        index_input += f" {cmd_arg} {tmp} "


with TemporaryDirectory() as tempdir:
    shell(
        "simpleaf index {extra} {index_input} "
        " --output {snakemake.output[0]} "
        " --threads {snakemake.threads} "
        " --work-dir {tempdir} "
        " {log} "
    )
