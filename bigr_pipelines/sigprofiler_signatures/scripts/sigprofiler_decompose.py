print('ok')

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

print('ok')

from SigProfilerExtractor import decomposition

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
        "output_dir", type=str, help="Path to output directory"
    )
    parser.add_argument(
        "-o", "--organism", type=str, help="Genome build",
        choices=["GRCh37", "GRCh38"], default="GRCh38"
    )
    parser.add_argument(
        "-v", "--verbose", help="Turn on verbose mode",
        action="store_true"
    )
    parser.add_argument(
        "--install", help="Install genome",
        action="store_true"
    )
    args = parser.parse_args()

    output_dir = f"{args.output_dir}/Res"
    logging.info(f"Working on {output_dir}, with {args.organism}")

    if args.install is True:
            genInstall.install('GRCh38')


    try:
        logging.info("Decomposing signatures for De novo SBS96")
        decomposition.decompose(
            f"{output_dir}/SBS96/Suggested_Solution/De_Novo_Solution/De_Novo_Solution_Signatures_SBS96.txt",
            f"{output_dir}/SBS96/Suggested_Solution/De_Novo_Solution/De_Novo_Solution_Activities_SBS96.txt",
            f"{output_dir}/SBS96/Samples.txt",
            f"{output_dir}/Deconvolution_SB96_DeNovo",
            genome_build=args.organism,
            verbose=args.verbose
        )
    except Exception as e:
        logging.error(e)
        raise

    try:
        logging.info("Decomposing signatures for Decomposed SBS96")
        decomposition.decompose(
            f"{output_dir}/SSBS96/Suggested_Solution/Decomposed_Solution/Decomposed_Solution_Signatures_SBS96.txt",
            f"{output_dir}/SBS96/Suggested_Solution/Decomposed_Solution/Decomposed_Solution_Activities_SBS96.txt",
            f"{output_dir}/SBS96/Samples.txt",
            f"{output_dir}/Deconvolution_SB96_Decomposed",
            genome_build=args.organism,
            verbose=args.verbose
        )
    except Exception as e:
        logging.error(e)
        raise

    try:
        logging.info("Decomposing signatures for De novo ID83")
        decomposition.decompose(
            f"{output_dir}/ID83/Suggested_Solution/De_Novo_Solution/De_Novo_Solution_Signatures_SBSINDEL.txt",
            f"{output_dir}/ID83/Suggested_Solution/De_Novo_Solution/De_Novo_Solution_Activities_SBSINDEL.txt",
            f"{output_dir}/ID83/Samples.txt",
            f"{output_dir}/Deconvolution_ID83_De_novo",
            genome_build=args.organism,
            verbose=args.verbose
        )
    except Exception as e:
        logging.error(e)
        raise

    try:
        logging.info("Decomposing signatures for Decomposed ID83")
        decomposition.decompose(
            f"{output_dir}/ID83/Suggested_Solution/Decomposed_Solution/Decomposed_Solution_Signatures_SBSINDEL.txt",
            f"{output_dir}/ID83/Suggested_Solution/Decomposed_Solution/Decomposed_Solution_Activities_SBSINDEL.txt",
            f"{output_dir}/ID83/Samples.txt",
            f"{output_dir}/Deconvolution_ID83_Decomposed",
            genome_build=args.organism,
            verbose=args.verbose
        )
    except Exception as e:
        logging.error(e)
        raise

    try:
        logging.info("Decomposing signatures for Decomposed DBS78")
        decomposition.decompose(
            f"{output_dir}/DBS78/Suggested_Solution/De_Novo_Solution/De_Novo_Solution_Signatures_SBSDINUC.txt",
            f"{output_dir}/DBS78/Suggested_Solution/De_Novo_Solution/De_Novo_Solution_Activities_SBSDINUC.txt",
            f"{output_dir}/DBS78/Samples.txt",
            f"{output_dir}/Deconvolution_DBS78_De_novo",
            genome_build=args.organism,
            verbose=args.verbose
        )
    except Exception as e:
        logging.error(e)
        raise


    try:
        logging.info("Decomposing signatures for Decomposed DBS78")
        decomposition.decompose(
            f"{output_dir}/DBS78/Suggested_Solution/Decomposed_Solution/Decomposed_Solution_Signatures_SBSDINUC.txt",
            f"{output_dir}/DBS78/Suggested_Solution/Decomposed_Solution/Decomposed_Solution_Activities_SBSDINUC.txt",
            f"{output_dir}/DBS78/Samples.txt",
            f"{output_dir}/Deconvolution_DBS78_Decomposed",
            genome_build=args.organism,
            verbose=args.verbose
        )
    except Exception as e:
        logging.error(e)
        raise
