#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""This wrapper estimates best solutions from SigProfilerExtractor results"""


import logging
import os

from os import path
from snakemake.utils import makedirs
from snakemake.shell import shell
from tempfile import TemporaryDirectory

from SigProfilerExtractor import estimate_best_solution as ebs

log = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)

logging.basicConfig(
    #filename=snakemake.log[0],
    filemode="w",
    level=logging.DEBUG
)

logging.info("Building expected directory structure")
with TemporaryDirectory() as tempdir:
    signatures_name = path.basename(snakemake.input["signatures"])

    os.symlink(
        snakemake.input["signatures"]
        path.join(tempdir, signatures_name)
    )

    ebs.estimate_solution(
        base_csvfile=snakemake.params.get("base_csvfile", "All_solutions_stat.csv"),
        All_solution=snakemake.params.get("all_solution", "All_Solutions"),
        genomes=snakemake.params.get("genomes", "Samples.txt"),
        output="results",
        **snakemake.params.get("estimate_extra", {})
    )

    logging.info("Retrieving results")
    csv = path.sep.join([
        tempfile,
        snakemake.params.get("all_solution", "All_Solutions"),
        "results",
        "All_solutions_stat.csv"
    ])

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
        "{csv} "
        "{snakemake.output['csv']} "
        "{log}"
    )

    pdf = path.sep.join([
        tempfile,
        snakemake.params.get("all_solution", "All_Solutions"),
        "results",
        "selection_plot.pdf"
    ])

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
        "{pdf} "
        "{snakemake.output['pdf']} "
        "{log}"
    )
