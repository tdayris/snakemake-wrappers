# coding: utf-8

"""Snakemake wrapper for macs2 callpeak"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2026, Thibault Dayris"
__license__ = "MIT"

from tempfile import TemporaryDirectory
from snakemake.shell import shell
from shlex import quote
from snakemake_wrapper_utils.snakemake import is_arg, move_files, list_arg, get_arg

log = snakemake.log_fmt_shell(stdout=True, stderr=True)
extra = snakemake.params.get("extra", "")

treatment = snakemake.input.get("treatment")
if treatment:
    if isinstance(treatment, list):
        treatment = [quote(str(t)) for t in treatment]
    else:
        treatment = quote(str(treatment))
    extra += f" --treatment {treatment}"

control = snakemake.input.get("control")
if control:
    if isinstance(control, list):
        control = [quote(str(c)) for c in control]
    extra += f" --control {control}"

barcodes = snakemake.input.get("barcodes")
if barcodes:
    extra += f" --barcodes {quote(barcodes)}"

with TemporaryDirectory() as tempdir:
    name = "NA"
    if is_arg("--name", extra):
        name = list_arg(extra)[get_arg("--name", extra) + 1]

    # Deal with required / optional output file(s)
    callpeak_results = f"{tempdir}/snake_outdir"
    expected_output_files = {
        "xls": f"{callpeak_results}/{name}_peaks.xls",
    }

    broad_peaks = snakemake.output.get("broad_peaks")
    gapped_peaks = snakemake.output.get("gapped_peaks")
    narrow_peaks = snakemake.output.get("narrow_peaks")
    if narrow_peaks and broad_peaks and gapped_peaks:
        raise ValueError("This tool can produce either narrow_peaks or broad_peaks")

    if narrow_peaks:
        # No option needed for narrow peak calls
        expected_output_files["narrow_peaks"] = str(
            f"{callpeak_results}/{name}_peaks.narrowPeak"
        )
        expected_output_files["summits"] = f"{callpeak_results}/{name}_summits.bed"
    elif broad_peaks or gapped_peaks:
        expected_output_files["broad_peaks"] = str(
            f"{callpeak_results}/{name}_peaks.broadPeak"
        )
        expected_output_files["gapped_peaks"] = str(
            f"{callpeak_results}/{name}_peaks.gappedPeak"
        )
        # Add option for broad peak calling
        extra += " --broad"

    control_bdg = snakemake.output.get("control_bdg")
    treatment_bdg = snakemake.output.get("treatment_bdg")
    if control_bdg or treatment_bdg:
        # control_bdg is produced if and only if control is provided
        if control:
            expected_output_files["control_bdg"] = str(
                f"{callpeak_results}/{name}_control_lambda.bdg"
            )
        expected_output_files["treatment_bdg"] = str(
            f"{callpeak_results}/{name}_treat_pileup.bdg"
        )
        # Add background option
        extra += " --bdg"

    # Move final file(s)
    move_cmds = f"tree -f {tempdir}  && "
    move_cmds += "&& ".join(move_files(snakemake, expected_output_files))

    shell(
        "( macs3 callpeak {extra} --tempdir {tempdir:q} "
        " --outdir {callpeak_results:q} && {move_cmds} ) {log}"
    )
