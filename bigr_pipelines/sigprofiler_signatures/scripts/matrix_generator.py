#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This script generates substitution matrices with SigProfiler

It expects all requirement (including namespace and directory architecture)
to be fullfiled.
"""

import argparse
import logging
import yaml


from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen

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
        "vcf", type=str, nargs=1, help="Path to input VCF file"
    )
    parser.add_argument(
        "yaml", type=str, nargs=1, help="Path to output yaml file"
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
    args = parser.parse_args()

    vcf_path = args.vcf
    out_yaml = args.yaml
    dir_name, vcf_name = vcf_path.rsplit("\t", 1)
    logging.info("Retrieving input information")
    logging.debug(f"Working at: {dir_name}")
    logging.debug(f"Working on: {vcf_name}")

    try:
        logging.info("Creating matrix")
        matrices = matGen.SigProfilerMatrixGeneratorFunc(
            project="test",
            genome=args.organism,
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

    logging.info(f"Saving results to {out_yaml}")
    with open(out_yaml, 'w') as yamlstream:
        yamlstream.write(yaml.dump(matrices, default_flow_style=False))
