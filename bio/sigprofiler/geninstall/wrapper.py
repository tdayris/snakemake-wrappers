#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""This wrapper runs SigProfiler over your VCF files"""


import logging
import os

from os import path
from tempfile import TemporaryDirectory

from SigProfilerMatrixGenerator import install as genInstall

logging.basicConfig(
    #filename=snakemake.log[0],
    filemode="w",
    level=logging.DEBUG
)

logging.info("Installing reference")
refname = (
    path.basename(snakemake.output["ref"][0][:-len(".gz")])
    if isinstance(snakemake.output["ref"], list)
    else path.basename(snakemake.output["ref"][:-len(".gz")])
)

if "ref" in snakemake.params.keys():
    genInstall.install(
        snakemake.params["ref"],
        rsync=snakemake.params.get("use_rsync", True),
        bash=not snakemake.params.get("use_rsync", True)
    )
elif "ref" in snakemake.input.keys():
    genInstall.install(
        refname,
        offline_files_path=snakemake.input["ref"]
    )
else:
    raise ValueError("Could not find any genome reference")
