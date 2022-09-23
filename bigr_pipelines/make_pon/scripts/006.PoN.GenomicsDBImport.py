#!/usr/bin/env python3
# coding: utf-8

import tempfile
from snakemake.shell import shell
from tempfile import TemporaryDirectory


log = snakemake.log_fmt_shell(stdout=False, stderr=True, append=True)

genomicsdb = snakemake.input.get("genomicsdb", "")
if genomicsdb:
    genomicsdb = f"--genomicsdb-update-workspace-path {genomicsdb}"
elif "genomicsdb" in snakemake.output.keys():
    genomicsdb = f"--genomicsdb-workspace-path {snakemake.output['genomicsdb']}"

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


variants = snakemake.input.get("variants")
if isinstance(variants, str):
    variants = f"--variant {variant}"
elif isinstance(variants, list):
    variants = "--variant " + "--variant".join(variants)
else:
    variants = ""

with tempfile.TemporaryDirectory() as tmpdir:
    shell(
        "gatk GenomicsDBImport "
        "--tmp-dir {tmpdir} "
        "--reference {snakemake.input.ref} "
        "{variants} "
        "{genomicsdb} "
        "{aln} "
        "{intervals} "
        "{intervals_list} "
        "{log} "
    )
