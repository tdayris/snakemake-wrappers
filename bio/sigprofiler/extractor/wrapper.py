#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""This wrapper rune SigProfiler Extractor over your vcf files"""

import logging
import os

from os import path
from snakemake.utils import makedirs
from snakemake.shell import shell
from tempfile import TemporaryDirectory

from SigProfilerExtractor import sigpro as sig

logging.basicConfig(
    #filename=snakemake.log[0],
    filemode="w",
    level=logging.DEBUG
)

input_type = None
if "vcf" in snakemake.input.keys():
    input_type = "vcf"
elif "matrix" in snakemake.input.keys():
    input_type = "matrix"
else:
    raise KeyError("Could not find either `vcf` or `matrix` in input keys")

logging.info("Building expected directory structure")
with TemporaryDirectory() as tempdir:
    input_data = None
    if "vcf" in snakemake.input.keys():
        original_vcf_dir = path.dirname(snakemake.input['vcf'][0])
        input_data = path.join(tempdir, path.basename(original_vcf_dir))
    else:
        input_data = snakemake.input["matrix"]

    project_name = snakemake.params.get("project_name", "Signatures")
    project_path = path.join(tempdir, project_name)

    makedirs(project_path)
    os.symlink(original_vcf_dir, input_data)

    logging.info("Extracting profiles")
    sig.sigProfilerExtractor(
        input_type=input_type,
        out_put=tempdir,
        input_data=input_data,
        reference_genome=snakemake.params.get("ref", "GRCh38"),
        opportunity_genome=snakemake.params.get("ref", "GRCh38"),
        cpu=snakemake.threads,
        gpu="gres" in snakemake.resources.keys(),
        **snakemake.params.get("extractor_extra", {})
    )

    logging.info("Retrieving results")
    shell(
        "rsync "
        "--verbose "
        "--checksum "
        "--recursive "
        "--update "
        "--links "
        "--partial "
        "--progress "
        "--human-readable "
        "{tempfile} "
        "{snakemake.output['signatures']} "
        "{log}"
    )
