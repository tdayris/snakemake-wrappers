#!/usr/bin/env python3
# coding: utf-8

import tempfile
from snakemake.shell import shell
from tempfile import TemporaryDirectory


log = snakemake.log_fmt_shell(stdout=False, stderr=True, append=True)


aln = snakemake.input.get("aln")
if isinstance(aln, str):
    aln = f"--input {aln}"
elif isinstance(aln, list):
    aln = "--input " + "--input ".join(aln)
else:
    aln = ""

intervals = snakemake.input.get("intervals")
if isinstance(intervals, str):
    intervals = f"--intervals {intervals}"
elif isinstance(intervals, list):
    intervals = "--intervals " + "--intervals".join(intervals)
else:
    intervals = ""

intervals_list = snakemake.output.get("interval_list", "")
if intervals_list:
    intervals_list = f"--output-interval-list-to-file {intervals_list}"



with tempfile.TemporaryDirectory() as tmpdir:
    shell(
        "gatk CreateSomaticPanelOfNormals "
        "--tmp-dir {tmpdir} "
        "--reference {snakemake.input.ref} "
        "--output {snakemake.output.pon} "
        "--variant {snakemake.input.gdb} "
        "{aln} "
        "{intervals} "
        "{log} "
    )
