__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2023, Thiault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

from snakemake.shell import shell
from distutils.spawn import find_executable
from os.path import dirname
from warnings import warn

import os

log = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)
extra = snakemake.params.get("extra", "")

# ALEVIN_FRY_HOME environment variable is required for
# simpleaf to run. It MUST be set.
alevin_fry_home = os.environ.get("ALEVIN_FRY_HOME")
if not alevin_fry_home:
    os.environ["ALEVIN_FRY_HOME"] = snakemake.params.get(
        "alevin_fry_home", dirname(snakemake.output[0])
    )
    warn(
        "ALEVIN_FRY_HOME variable not found in environment. "
        f"It has been set to {os.environ['ALEVIN_FRY_HOME']}"
    )


# Searching for executables if user did not
# provide a path as input.
salmon = snakemake.input.get("salmon", find_executable("salmon"))
piscem = snakemake.input.get("piscem", find_executable("piscem"))
alevin = snakemake.input.get("alevin_fry", find_executable("alevin-fry"))
pyroe = snakemake.input.get("pyroe", find_executable("pyroe"))


shell(
    "simpleaf set-paths {extra} "
    "--salmon {salmon} "
    "--piscem {piscem} "
    "--alevin-fry {alevin} "
    "--pyroe {pyroe} "
    "{log}"
)

if snakemake.output[0] != "simpleaf_info.json":
    shell("mv --verbose simpleaf_info.json {snakemake.output[0]} {log}")