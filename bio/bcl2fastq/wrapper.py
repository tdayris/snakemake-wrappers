#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""This is the Snakemake Wrapper for bcl2fastq"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2020, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"


from snakemake.shell import shell
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

# extra parameters
extra = snakemake.params.get("extra", "")

## Input parameters
# path to the sample sheet
sample_sheet = ""
if "sample_sheet" in snakemake.input.keys():
    sample_sheet = "--sample-sheet {}".format(snakemake.input["sample_sheet"])

# path to input directory
input_dir = ""
if "input_dir" in snakemake.input.keys():
    input_dir = "--input-dir {}".format(snakemake.input["input_dir"])

# path to runfolder directory
run_dir = ""
if "run_dir" in snakemake.params.keys():
    run_dir = "--runfolder-dir {}".format(snakemake.input["run_dir"])

# path to intensities directory, requires an input dir (see above)
intensities = ""
if "intensities" in snakemake.input.keys():
    if "input_dir" not in snakemake.input.keys():
        raise KeyError(
            "If intensities is specified, input_dir must also be specified."
        )
    intensities = "--intensities-dir {}".format(snakemake.input["intensities"])


## Output parameters
# path to demultiplexed output
out_dir = ""
if "out_dir" in snakemake.output.keys():
    out_dir = " --output-dir {}".format(snakemake.output["out_dir"])

# path to demultiplexing statistics directory
interop_dir = ""
if "interop_dir" in snakemake.output.keys():
    interop_dir = "--interop-dir {}".format(snakemake.output["interop_dir"])

# path to human-readable demultiplexing statistics directory
stats = ""
if "stats_dir" in snakemake.output.keys():
    stats = "--stats-dir {}".format(snakemake.output["stats_dir"])

# path to reporting directory
reports = ""
if "reports_dir" in snakemake.output.keys():
    reports = "--reports-dir {}".format(snakemake.output["reports_dir"])


# The number of threads shall be spreaded between loading, processing and
# writing. Knowing that the processing step requires more threads than the
# other ones. Documentation advises on twice more threads while on the
# processing step.
max_threads = snakemake.threads

if max_threads < 3:
    raise ValueError("At least four threads are required for bcl2fastq.")
elif max_threads == 3:
    # In case of only three threads requested, the code below would fail
    # by returning 0 threds in loading/writing.
    reserved_threads = [1, 1, 1]
else:
    # Optimal division according to documentation, knowing that all number
    # must be plain integers.
    reserved_threads = [
        int(max_threads/4),
        int(max_threads/2),
        int(max_threads/4)
    ]

# With odd number, the above will result in less reserved threads than
# the number provided by the user. Below, we simply add the difference
# (if any) to the processing step
reserved_threads[1] += max_threads - sum(reserved_threads)

shell(
    "bcl_to_fastq "
    " {input_dir} "  # path to input directory
    " {intensities} "  # path to intensities directory
    " {run_dir} "  # path to runfolder directory
    " {sample_sheet} "  # path to the sample sheet
    " {out_dir} "  # path to demultiplexed output
    " {interop_dir} "  # path to demultiplexing statistics directory
    " {stats_dir} "  # path to human-readable statistics directory
    " {reports_dir} "  # path to reporting directory
    # number of threads used for loading BCL data
    " --loading-threads {reserved_threads[0]} "
    # number of threads used for processing demultiplexed data
    " --processing-threads {reserved_threads[1]} "
    # number of threads used for writing FASTQ data
    " --writing-threads {reserved_threads[2]} "
    " {extra} "  # Extra parameters
    " {log} "  # Logging
)
