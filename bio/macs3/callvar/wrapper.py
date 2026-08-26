# coding: utf-8

"""Snakemake wrapper for macs2 callvar"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2026, Thibault Dayris"
__license__ = "MIT"

from snakemake.shell import shell
from snakemake.ioutils import subpath
from shlex import quote
from snakemake_wrapper_utils.bcftools import get_bcftools_opts, infer_out_format

log = snakemake.log_fmt_shell(stdout=True, stderr=True)
extra = snakemake.params.get("extra", "")

# Optional input file path
control = snakemake.input.get("control")
if control:
    extra += f" --control {quote(control)}"

# Handle output file path
outfmt = infer_out_format(str(snakemake.output.call))
out_cmd = ""
if outfmt == "v":
    outdir = quote(subpath(str(snakemake.output.call), parent=True))
    outfile = quote(subpath(str(snakemake.output.call), basename=True))
    extra += f" --outdir {outdir} --ofile {outfile}"
else:
    extra += f" --ofile '/dev/stdout'"
    bcftools_opts = get_bcftools_opts(snakemake, parse_memory=False)
    out_cmd = f"| bcftools view {bcftools_opts}"

call_idx = snakemake.output.get("call_idx")
if call_idx:
    index_extra=""
    if str(call_idx).endswith("csi"):
        index_extra += " --csi"
    elif str(call_idx).endswith("tbi"):
        index_extra += " --tbi"
    out_cmd += f" && bcftools index {index_extra} {quote(str(snakemake.output.call))}"

shell(
    "( macs3 callvar {extra} --peak {snakemake.input.peak:q} "
    "--multiple-processing {snakemake.threads} "
    "--treatment {snakemake.input.treatment:q} {out_cmd} ) {log} "
)
