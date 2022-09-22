#!/usr/bin/env python3
# coding: utf-8

from typing import Dict, List

def get_t2n_extra(database: str) -> str:
    """Get the correct awk print command to extract human readable term 2 name"""
    two_cols = [
        "WikiPathway",
        "Reactome",
        "Human_Phenotype_Onthology",
        "Human_Protein_Atlas",
    ]
    first_col_readable = [
        "MSigDB_c1_Positional_Gene_Sets",
        "MSigDB_c2_Curated_Gene_Sets",
        "MSigDB_c3_Regulatory_Targets",
        "MSigDB_c4_Cancer_Oriented_MicroArray",
        "MSigDB_c5_Onthology_Gene_Sets",
        "MSigDB_c6_Oncogenic_Signatures",
        "MSigDB_c7_Immunologic_Signatures",
        "MSigDB_c8_Cell_Type_Signatures",
        "MSigDB_hallmark",
        "Chemical_genetic_perturbation",
        "Canonical_pathways",
        "Biocarta",
        "KEGG",
        "PID",
        "MSigDB_Reactome",
        "MSigDB_WikiPathways",
        "MicroRNA_targets",
        "miRDB",
        "MicroRNA_legacy",
        "Transcription_factor_targets",
        "Transcription_factor_GTRD_Kolmykov",
        "Transcription_factor_targets_legacy",
        "Cancer_gene_neighborhood",
        "Cancer_modules",
        "Gene_Othology_all",
        "Gene_Othology_BP",
        "Gene_Othology_CC",
        "Gene_Othology_MF",
        "MSigDB_Human_Phenotype_Onthology",
        "ImmuneSigDB",
        "Vaccine_Response",
        "WikiPathway",
        "MiRNABase",
        "CORUM",
    ]

    if database in two_cols:
        return '{print $1"\t"$2}'
    
    if database in first_col_readable:
        return '{print $1"\t"$1}'

    return '{print $1"\t"$1}'


def db_keytype(gmts: Dict[str, str], ppis: Dict[str, str]) -> Dict[str, str]:
    """get correct keytype for a given database"""
    result = {}
    for db_name, db_path in gmts.items():
        if db_path.lower().endswith(".entrez.gmt"):
            result[db_name] = "ENTREZID"
        elif db_path.lower().endswith(".ensg.gmt"):
            result[db_name] = "ENSEMBL"
        elif db_path.lower().endswith("symbols.gmt"):
            result[db_name] = "SYMBOL"
        else:
            result[db_name] = "ENTREZID"

    for ppi in ppis.keys():
        result[ppi] = "ENSEMBLPROT"

    return result


def plot_list(plots: List[str], 
              methods: List[str],
              comparisons: List[str], 
              db_key: Dict[str, str]) -> List[str]:
    """Return list of plots with correct key type and database"""
    result = [
        f"data_output/{comparison}/{db_name}.{keytype}/{plot}.{method}.png"
        for comparison in comparisons
        for db_name, keytype in db_key.items()
        for plot in plots
        for method in methods
    ]
    result += [
        f"data_output/all_comparisons/{db_name}.{keytype}.png"
        for db_name, keytype in db_key.items()
    ]
    return result

def tsv_list(db_key: Dict[str, str], comparisons: List[str]) -> List[str]:
    """Return list of expected TSV files"""
    return [
        f"results/{comparison}/{db_name}.{keytype}/enrich.{comparison}.{keytype}.tsv"
        for comparison in comparisons
        for db_name, keytype in db_key.items()
    ]