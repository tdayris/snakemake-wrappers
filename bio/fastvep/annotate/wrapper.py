# coding: utf-8

"""Snakemake wrapper for fastvep external annotation source build"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2026, Thibault Dayris"
__license__ = "MIT"

from os.path import commonprefix
from snakemake.shell import shell
from snakemake_wrapper_utils.snakemake import get_format, is_arg
from shlex import quote

extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

out_fmt = get_format(snakemake.output.variants)
if out_fmt in ("vcf", "json"):
    extra += f" --output-format {quote(out_fmt)}"
elif out_fmt == "tsv":
    extra += f" --output-format 'tab'"

gff = snakemake.input.get("gff")
if gff:
    extra += f" --gff3 {quote(gff)}"

fasta = snakemake.input.get("fasta")
if fasta:
    extra += f" --fasta {quote(fasta)}"

cache = snakemake.input.get("cache")
if cache:
    extra += f" --cache-dir {quote(cache)}"

transcript_cache = snakemake.input.get("transcript_cache")
if transcript_cache:
    extra += f" --transcript-cache {quote(transcript_cache)}"

sa_dir = snakemake.input.get("sa_dir")
if sa_dir:
    extra += f" --sa-dir {quote(sa_dir)}"

acmg_config = snakemake.input.get("acmg_config")
if acmg_config:
    extra += f" --acmg-config {quote(acmg_config)}"
    if not is_arg("--acmg", extra):
        extra += " --acmg"

functional_evidence = snakemake.input.get("functional_evidence")
if functional_evidence:
    extra += f" --functional-evidence {quote(functional_evidence)}"

panel = snakemake.input.get("panel")
if panel:
    extra += f" --gene-list {quote(panel)}"

qc_rules = snakemake.input.get("qc_rules")
if qc_rules:
    extra += f" --qc-rules {quote(qc_rules)}"


shell(
    "fastvep annotate {extra} "
    "--input {snakemake.input.vcf:q} "
    "--output {snakemake.output.vcf:q} "
    "{log} "
)
