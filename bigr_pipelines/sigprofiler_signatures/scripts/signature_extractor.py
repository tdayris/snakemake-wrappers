#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This script extracts signatures fom matrices

It expects all requirement (including namespace and directory architecture)
to be fullfiled.
"""

import argparse
import logging
import sys

from SigProfilerMatrixGenerator import install as genInstall
from SigProfilerExtractor import sigpro as sig
#from SigProfilerExtractor import estimate_best_solution as ebs


logging.basicConfig(
    #filename=snakemake.log[0],
    filemode="w",
    level=logging.DEBUG
)


if __name__ == '__main__':
    logging.info("Parsing command line input")
    parser = argparse.ArgumentParser(
        description=sys.modules[__name__].__doc__,
        epilog="This script does not perform any magic. Check the results."
    )

    parser.add_argument(
        "input", type=str, help="Path to input directory"
    )
    parser.add_argument(
        "output", type=str, help="Path to output directory"
    )
    parser.add_argument(
        "-t", "--threads", type=int, help="Maximum number of threads",
        default=1
    )
    parser.add_argument(
        "-o",
        "--organism",
        type=str,
        help="Organism identifier",
        choices=["GRCh38", "GRCh37"],
        default="GRCh38"
    )
    parser.add_argument(
        "--install", help="Install genome",
        action="store_true"
    )
    parser.add_argument(
        "-g",
        "--gpu",
        default=False,
        action="store_true"
    )
    args = parser.parse_args()

    input_dir = args.input
    out_dir = args.output
    logging.info("Retrieving command line information")
    logging.debug(
        f"Working on {input_dir}, ({args.organism})"
        f"with {args.threads} thread(s), "
        f"GPU:{args.gpu}, "
        f"Saving results in {out_dir}"
    )

    #try:
    #    if args.install is True:
    #        genInstall.install('GRCh38')

    #    logging.info("Extracting signatures")
    #    sig.sigProfilerExtractor(
    #        input_type="vcf",
    #        output=out_dir,
    #        input_data=input_dir,
    #        reference_genome=args.organism,
    #        cpu=args.threads,
    #        gpu=args.gpu,
            #min_stability=0.1
    #    )

    #    ebs.estimate_solution(base_csvfile="All_solutions_stat.csv",
    #      All_solution="All_Solutions",
    #      genomes="Samples.txt",
    #      output="results",
    #      title="Selection_Plot",
    #      stability=0.8,
    #      min_stability=0.2,
    #      combined_stability=1.25)
    #except Exception as e:
    #    logging.error(e)
    #    raise
