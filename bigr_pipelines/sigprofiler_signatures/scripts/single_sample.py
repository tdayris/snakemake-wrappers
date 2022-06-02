#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This script extracts signatures fom matrices

It expects all requirement (including namespace and directory architecture)
to be fullfiled.
"""

import argparse
import logging
import os.path
import pandas
import sys

from SigProfilerMatrixGenerator import install as genInstall
from sigproSS import spss


logging.basicConfig(
    #filename=snakemake.log[0],
    filemode="w",
    level=logging.DEBUG
)


if __name__ == '__main__':
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
        "--sbs", type=str, help="Path to SBS CSV file", default=None
    )

    parser.add_argument(
        "--dbs", type=str, help="Path to DBS CSV file", default=None
    )

    parser.add_argument(
        "--id", type=str, help="Path to ID CSV file", default=None
    )

    # parser.add_argument(
    #     "--tbs",  type=str, help="Path to TSB CSV file", default=None
    # )

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

    args = parser.parse_args()

    logging.basicConfig(
        #filename=args.logging,
        filemode="w",
        level=logging.DEBUG
    )

    sig_database = "default"
    check_rules = True
    if args.sbs is not None:
        sig_database = pandas.read_csv(
            args.dbs, sep=",", header=0, index_col=None
        )
        sig_database["Mutation Type"] = [
            subtype[0] + "[" + sub + "]" + subtype[1]
            for subtype, sub in zip(sig_database["SubType"], sig_database["Type"])
        ]
        sig_database.set_index("Mutation Type", inplace=True)
        del sig_database["SubType"]
        del sig_database["Type"]
        logging.info("Using SBS signature profiles")
        logging.debug(sig_database.head())
    if args.dbs is not None:
        sig_database = pandas.read_csv(
            args.dbs, sep=",", header=0, index_col=0
        )
        logging.info("Using DBS signature profiles")
        logging.debug(sig_database.head())
        check_rules = False
    elif args.id is not None:
        sig_database = pandas.read_csv(
            args.id, sep=",", header=0, index_col=None
        )
        sig_database["Mutation Type"] = [
            i.replace("+", "p") for i in sig_database["Mutation Type"]
        ]
        sig_database.set_index("Mutation Type", inplace=True)
        logging.info("Using ID signature profiles")
        logging.debug(sig_database.head())
        check_rules = False
    else:
        logging.info("Using default SBS signature profiles")


    input_dir = os.path.dirname(args.input)
    out_dir = args.output
    logging.info("Retrieving command line information")
    logging.debug(
        f"Working on {args.input} (moved in input_dir) "
        f"Saving results in {out_dir} (analysis for {args.organism})"
    )

    try:
        if args.install is True:
            logging.info("Installing genome")
            logging.warning(
                "If multiple instences of SigProfile simultaneously install"
                " genome, the result will either fail, or end corrupted."
            )
            genInstall.install('GRCh38')

        logging.info(f"Analysing {args.input}")
        spss.single_sample(
            input_dir,
            out_dir,
            ref="GRCh38",
            sig_database=sig_database,
            check_rules=check_rules,
            exome=True
         )
    except Exception as e:
        logging.error(e)
        raise
