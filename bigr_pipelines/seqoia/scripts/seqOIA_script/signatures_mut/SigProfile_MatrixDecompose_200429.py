#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 29 09:29:13 2020

@author: mneou
"""
if __name__ == "__main__":
    import sys
    dossier_output = sys.argv[1]
    print(dossier_output)
    ## Extraction des signatures  ## 
    print('extraction2')
    from SigProfilerExtractor import decomposition as decomp
    ##### SBS96 #####
    signatures = dossier_output+ "SBS96/Suggested_Solution/De_Novo_Solution/De_Novo_Solution_Signatures_SBS96.txt"
    activities=dossier_output + "SBS96/Suggested_Solution/De_Novo_Solution/De_Novo_Solution_Activities_SBS96.txt"
    samples=dossier_output + "SBS96/Samples.txt"
    output= dossier_output + "Deconvolution_SB96_DeNovo"
    decomp.decompose(signatures, activities, samples, output, genome_build="GRCh38", verbose=True)
    ####
    signatures = dossier_output+ "SBS96/Suggested_Solution/Decomposed_Solution/Decomposed_Solution_Signatures_SBS96.txt"
    activities=dossier_output + "SBS96/Suggested_Solution/Decomposed_Solution/Decomposed_Solution_Activities_SBS96.txt"
    #####
    samples=dossier_output + "SBS96/Samples.txt"
    output= dossier_output + "Deconvolution_SB96_Decomposed"
    decomp.decompose(signatures, activities, samples, output, genome_build="GRCh38", verbose=False)
    ##### ID83 #####
    signatures = dossier_output+ "ID83/Suggested_Solution/Decomposed_Solution/Decomposed_Solution_Signatures_SBSINDEL.txt"
    activities=dossier_output + "ID83/Suggested_Solution/Decomposed_Solution/Decomposed_Solution_Activities_SBSINDEL.txt"
    samples=dossier_output + "ID83/Samples.txt"
    output= dossier_output + "Deconvolution_ID83_Decomposed"
    decomp.decompose(signatures, activities, samples, output, genome_build="GRCh38", verbose=False)
    ###
    signatures = dossier_output+ "ID83/Suggested_Solution/De_Novo_Solution/De_Novo_Solution_Signatures_SBSINDEL.txt"
    activities=dossier_output + "ID83/Suggested_Solution/De_Novo_Solution/De_Novo_Solution_Activities_SBSINDEL.txt"
    samples=dossier_output + "ID83/Samples.txt"
    output= dossier_output + "Deconvolution_ID83_De_novo"
    decomp.decompose(signatures, activities, samples, output, genome_build="GRCh38", verbose=False)
    ##### DBS78 #####   
    signatures = dossier_output+ "DBS78/Suggested_Solution/Decomposed_Solution/Decomposed_Solution_Signatures_SBSDINUC.txt"
    activities=dossier_output + "DBS78/Suggested_Solution/Decomposed_Solution/Decomposed_Solution_Activities_SBSDINUC.txt"
    samples=dossier_output + "DBS78/Samples.txt"
    output= dossier_output + "Deconvolution_DBS78_Decomposed"
    decomp.decompose(signatures, activities, samples, output, genome_build="GRCh38", verbose=False)
    ###
    signatures = dossier_output+ "DBS78/Suggested_Solution/De_Novo_Solution/De_Novo_Solution_Signatures_SBSDINUC.txt"
    activities=dossier_output + "DBS78/Suggested_Solution/De_Novo_Solution/De_Novo_Solution_Signatures_SBSDINUC.txt"
    samples=dossier_output + "DBS78/Samples.txt"
    output= dossier_output + "Deconvolution_DBS78_De_novo"
    decomp.decompose(signatures, activities, samples, output, genome_build="GRCh38", verbose=False)
