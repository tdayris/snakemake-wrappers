# coding: utf-8

"""Snakemake wrapper for running quality control with fastcat qc."""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2026, Thibault Dayris"
__license__ = "MIT"

from os.path import isdir
from shlex import quote
from snakemake.shell import shell
from snakemake_wrapper_utils.snakemake import get_format, is_arg, move_files
from tempfile import TemporaryDirectory

extra = snakemake.params.get("extra", "")
crabz_extra = snakemake.params.get("crabz_extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)
expected_output_files = {}

# Input can be either a single/multiple fastq, or a directory of fastq files,
# or a single/multiple BAM, or a directory of BAM files, 
# or a single/multiple file listing input files, or a directory of files listing inputs
input_formats = []
for file_path in snakemake.input.reads:
    # Avoid typing issues
    file_path = str(file_path)
    input_formats.append(get_format(file_path))
    if isdir(file_path) and not is_arg("--recurse", extra):
        extra += " --recurse"

# Detecting input format if possible
input_formats = list(set(input_formats))
if len(input_formats) > 1:
    raise ValueError("Provide either fastqs, bams, OR txt. Not a mix of them.")
elif input_formats[0] == "fastq":
    extra += f" --input-fmt='fastq'"
elif input_formats[0] == "bam":
    extra += f" --input-fmt='bam'"

# Optional IO file(s)
bed = snakemake.input.get("bed")
if bed:
    extra += f" --bed {quote(bed)}"

json_progress = snakemake.output.get("json_progress")
if json_progress:
    extra += f" --jsonl-progress {quote(json_progress)}"

coverage_bed = snakemake.input.get("coverage_bed")
if coverage_bed:
    extra += f" --coverage-beds {coverage_bed} "

reads = snakemake.output.get("reads")
out_cmd = ""
# Fastq files are concatenated and sent to standard output
if reads:
    if input_formats[0] == "bam":
        raise ValueError("No reads output produced with bam input files.")
    out_cmd = f" > {quote(str(reads))}"
    if str(reads).endswith(".gz"):
        out_cmd = f" | crabz > {quote(str(reads))}"


with TemporaryDirectory() as tempdir:
    demux = snakemake.output.get("demux")
    if demux:
        expected_output_files["demux"] = f"{tempdir}/barcode*"
        extra += f" --demultiplex"
    coverage_all = snakemake.output.get("coverage")
    if coverage:
        expected_output_files["coverage"] = f"{tempdir}/coverage"
        extra += f" --coverage"

    per_read = snakemake.output.get("per_read")
    if per_read:
        extra += f" --per-read"
        expected_output_files["per_read"] = f"{tempdir}/per_read.tsv"

    html = snakemake.output.get("html")
    if html:
        extra += f" --html-report"
        expected_output_files["html"] = f"{tempdir}/fastcat_report.html"

    move_cmds = " && ".join(move_files(snakemake, expected_output_files))

    shell(
        "( fastcat qc --threads {snakemake.threads} {extra} {snakemake.input.reads} "
        "  --output '{tempdir}/fastcat_qc' {out_cmd} && {move_cmds} ) {log}"
    )



