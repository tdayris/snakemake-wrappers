#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Annotate VCF file with VCFtools perl script"""


from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=False, stderr=True)
extra = snakemake.params.get("extra", "")

if "description" in snakemake.input.keys():
    extra += " --description {}".format(snakemake.input["description"])

if "annotation" in snakemake.input.keys():
    extra += " --annotations {}".format(snakemake.input["annotation"])

stream = "gunzip -c" if str(snakemake.input.vcf).endswith(".gz") else "cat"

sort = "vcf-sort | " if snakemake.params.get("sort", True) is True else ""

if str(snakemake.output.vcf).endswith(".gz"):
    if snakemake.threads != 4:
        raise ValueError("This wrapper requires exactly three threads to work")
    shell(
        "({stream} {snakemake.input.vcf} | "
        "vcf-annotate {extra} | "
        "{sort} "
        "pbgzip -c) > {snakemake.output.vcf} {log}"
    )
else:
    if snakemake.threads != 3:
        raise ValueError("This wrapper requires exactly two threads to work")
    shell(
        "({stream} {snakemake.input.vcf} | "
        "{sort} "
        "vcf-annotate {extra}) "
        "> {snakemake.output.vcf} {log}"
    )