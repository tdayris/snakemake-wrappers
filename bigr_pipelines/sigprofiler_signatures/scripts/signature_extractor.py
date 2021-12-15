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

from SigProfilerExtractor import sigpro as sig


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
        "input", type=str, nargs=1, help="Path to input directory"
    )
    parser.add_argument(
        "-t", "--threads", type=int, nargs=1, help="Maximum number of threads",
        default=1
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
        "-g",
        "--gpu",
        default=False,
        action="store_true"
    )
    args = parser.parse_args()

    input_dir = args.input
    out_dir = f"{input_dir}/Res"
    logging.info("Retrieving command line information")
    logging.debug(
        f"Working on {input_dir}, ({args.organism})"
        f"with {args.threads} thread(s), "
        f"GPU:{args.gpu}"
    )

    try:
        logging.info("Extracting signatures")
        sig.sigProfilerExtractor(
            "vcf",
            dossier_output,
            dossier_input,
            args.organism,
            cpu=args.threads,
            gpu=args.gpu
        )
    except Exception as e:
        logging.error(e)
        raise
