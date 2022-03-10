#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This script generates substitution matrices with SigProfiler

It expects all requirement (including namespace and directory architecture)
to be fullfiled.
"""

import argparse
import logging
import sys

from SigProfilerMatrixGenerator import install as genInstall
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen

logging.basicConfig(
    #filename=snakemake.log[0],
    filemode="w",
    level=logging.DEBUG
)

logging.getLogger('matplotlib.font_manager').disabled = True
logging.getLogger('matplotlib.backends.backend_pdf').disabled = True

if __name__ == '__main__':
    logging.info("Parsing command line input")
    parser = argparse.ArgumentParser(
        description=sys.modules[__name__].__doc__,
        epilog="This script does not perform any magic. Check the results."
    )

    parser.add_argument(
        "vcf", type=str, help="Path to input VCF file"
    )
    parser.add_argument(
        "-o",
        "--organism",
        type=str,
        nargs=1,
        help="Organism identifier",
        choices=["GRCh38", "GRCh37"],
        default="GRCh38"
    )
    parser.add_argument(
        "-s", "--sample-name",
        help="Name of the sample to process",
        type=str,
        default=None
    )
    parser.add_argument(
        "--install", help="Install genome",
        action="store_true"
    )
    args = parser.parse_args()

    vcf_path = args.vcf
    dir_name, vcf_name = vcf_path.rsplit("/", 1)
    organism = args.organism[0] if isinstance(args.organism, list) else args.organism
    organism = str(organism)
    sample_name = args.sample_name if args.sample_name is not None else vcf_name.split(".")[0]
    logging.info("Retrieving input information")
    logging.debug(f"Working at: {dir_name}")
    logging.debug(f"Working on: {vcf_name}")
    
    try:
        if args.install is True:
            genInstall.install('GRCh38')

        logging.info("Creating matrix")
        matrices = matGen.SigProfilerMatrixGeneratorFunc(
            project="test",
            genome=organism,
            vcfFiles=dir_name,
            exome=True,
            bed_file=None,
            chrom_based=None,
            plot=True,
            tsb_stat=False,
            seqInfo=False,
            cushion=1000
        )
    except Exception as e:
        logging.error(e)
        raise

