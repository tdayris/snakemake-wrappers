#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""This wrapper runs SigProfiler MatrixGenerator over your VCF files"""


import logging
import os

from os import path
from snakemake.utils import makedirs
from snakemake.shell import shell
from tempfile import TemporaryDirectory

from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen

log = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)

logging.basicConfig(
    #filename=snakemake.log[0],
    filemode="d",
    level=logging.DEBUG
)

def list_files(startpath):
    for root, dirs, files in os.walk(startpath):
        level = root.replace(startpath, '').count(os.sep)
        indent = ' ' * 4 * (level)
        logging.debug('{}{}/'.format(indent, os.path.basename(root)))
        subindent = ' ' * 4 * (level + 1)
        for f in files:
            logging.debug('{}{}'.format(subindent, f))

logging.info("Building expected directory structure")
with TemporaryDirectory() as tempdir:
    vcf_basename = path.basename(snakemake.input["vcf"])
    vcf_path = path.join(tempdir, vcf_basename)

    project_name = snakemake.params.get("project_name", "Signatures")
    project_path = path.join(tempdir, project_name)

    refname = path.basename(snakemake.input["ref"][:-len(".tar.gz")])

    makedirs(project_path)
    os.symlink(
        snakemake.input["vcf"],
        vcf_path
    )
    list_files(tempdir)


    logging.info("Building matrices")
    matrices = matGen.SigProfilerMatrixGeneratorFunc(
        project_name,
        refname,
        tempdir,
        **snakemake.params.get(
            "matrix_extra",
            {
                "plot": True,
                "exome": False,
                "bed_file": None,
                "chrom_based": False,
                "tsb_stat": False,
                "seqInfo": False,
                "cushion": 100
            }
        )
    )
    logging.debug(matrices)

    logging.info("Retrieving results")
    list_files(tempdir)

    shell(
        "rsync "
        "--verbose "
        "--checksum "
        "--recursive "
        "--update "
        "--links "
        "--human-readable "
        "{project_path} "
        "{snakemake.output['matrices']} "
        "{log}"
    )
    list_files(snakemake.output['matrices'])
