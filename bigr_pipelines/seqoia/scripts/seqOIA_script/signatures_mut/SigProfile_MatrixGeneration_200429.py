#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 29 09:27:53 2020

@author: mneou
"""

if __name__ == "__main__":
    from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
    import argparse
    import sys
    import os
    parser = argparse.ArgumentParser(description='Mutational signature matrix generation')
    parser.add_argument('-i', '--input', type=str, nargs=1,action='append', help='Input files path', required=True)
    args = parser.parse_args()
    vcf_input = args.input[0][0]
    #########
    print("input :", vcf_input)
    print("folder :", os.path.dirname(vcf_input))

    matrices = matGen.SigProfilerMatrixGeneratorFunc("test", "GRCh38", os.path.dirname(vcf_input), plot=True, exome=True, bed_file=None, chrom_based=None, tsb_stat=False, seqInfo=False, cushion=1000)
    ##
   
