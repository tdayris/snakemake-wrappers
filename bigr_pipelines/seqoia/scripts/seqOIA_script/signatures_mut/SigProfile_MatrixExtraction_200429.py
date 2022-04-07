#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 29 09:29:13 2020

@author: mneou
"""
if __name__ == "__main__":

    
    import sys    
    dossier_input = '/'.join(sys.argv[1].split("/")[0:-1])

    
    dossier_output =  dossier_input + "/Res"
    print(dossier_input)
    print(dossier_output)
    ## Extraction des signatures  ## 
    print('extraction2')
    from SigProfilerExtractor import sigpro as sig
    #dossier_output = "VCFs/Results"
    #dossier_input = "VCFs"
    # 
    sig.sigProfilerExtractor("vcf", dossier_output, dossier_input,"GRCh38", cpu=8)
    
    print("fin")
