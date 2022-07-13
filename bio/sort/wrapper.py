#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from snakemake.shell import shell


log = snakemake.log_fmt_shell(stdout=False, stderr=True)
tmpdir = snakemake.resources.get("tmpdir", "/tmp")

extra = snakemake.params.get("extra", "")

keys = snakemake.params.get("columns", [])
if isinstance(keys, str) or isinstance(keys, int):
    keys = f"-k{keys}"
elif isinstance(keys, list):
    keys = " ".join([f"-k{k}" for k in keys])
else:
    keys = ""

shell(
    "sort --parallel={snakemake.threads} "
    "--temporary-directory={tmpdir} "
    "{keys} {extra} {snakemake.input[0]} "
    "> {snakemake.output[0]} {log}"
)
