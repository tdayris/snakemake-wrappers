#!/usr/bin/python3.11
# -*- coding: utf-8 -*-

import re

from typing import Optional, Tuple

def get_snpeff_ann(info: str) -> Optional[str]:
    """From a VCF info field, return SnpEff ANN value"""
    regex = r'ANN=([^;]+)'
    try:
        return re.match(regex, info)[1]
    except IndexError:
        return None


def get_variant_class(info: str) -> str:
    """From a VCF info field, return SnpSift Variant Class"""
    if "SNP;" in info:
        return "SNP"
    if "DEL;" in info:
        return "DEL"
    if "INS" in info:
        return "INS"
    if "MNP" in info:
        return "MNP"
    if "MIXED" in info:
        return "MIXED"
    return ""


def get_vep_ranking(
    rank: Optional[str] = None,
    annotation: Optional[str] = None
) -> Tuple[str]:
    """From a ranking and an annotation, return the correct vep information"""
    if "intron" in annotation:
        return "", rank
    return rank, ""


def get_pos(ratio: Optional[str] = None) -> Tuple[str]:
    """From a cdna_pos/cdna_len, return only the first"""
    if "/" in ratio:
        return ratio.split("/")[0]
    return ""


def getaminoacid(hgvsp: Optional[str] = None)


def snpeff2vep(ann: str, rs: str, variant_class: str) -> str:
    """Convert SnpEff ANNotation field to a VEP one"""
    ann = ann.split("|")
    ann_dict = {
        "Allele": ann[0],
        "Annotation": ann[1],
        "Putative_impact": ann[2],
        "Gene Name": ann[3],
        "Gene ID": ann[4],
        "Feature type": ann[5],
        "Feature ID": ann[6],
        "Transcript biotype": ann[7],
        "Rank/Total": ann[8],
        "HGVS.c": ann[9],
        "HGVS.p": ann[10],
        "cDNA_position/cDNA_len": ann[11],
        "CDS_position/CDS_len": ann[12],
        "Protein_position/Protein_len": ann[13],
        "Distance": ann[14],
        "Comments": ann[15],
    }

    exon, intron = get_vep_ranking(
        rank=ann_dict["Rank/Total"].split("/"),
        annotation=ann_dict["Annotation"],
    )
    cdna = get_pos(ratio=ann_dict["cDNA_position/cDNA_len"],)
    cds = get_pos(ratio=ann_dict["CDS_position/CDS_len"],)
    protein = get_pos(ratio=ann_dict["Protein_position/Protein_len"],)
    aa = get_amino_acid(hgvsp=ann_dict["HGVS.p"], hgvsc=ann_dict["HGVS.c"])

    vep_dict = {
        "Allele": ann_dict["Allele"],
        "Consequence": ann_dict["Annotation"],
        "IMPACT": ann_dict["Putative_impact"],
        "SYMBOL": ann_dict["Gene Name"],
        "Gene": ann_dict["Gene ID"],
        "Feature_type": ann_dict["Feature type"],
        "Feature": ann_dict["Feature ID"],
        "BIOTYPE": ann_dict["Transcript biotype"],
        "EXON": exon,
        "INTRON": intron,
        "HGVSc": ann_dict["HGVS.c"],
        "HGVSp": ann_dict["HGVS.p"],
        "cDNA_position": cdna,
        "CDS_position": cds,
        "Protein_position": protein,
        "Amino_acids": "", # FILL? Should be like G/D (AminoAcidOneLetterReference/AminoAcidOneLetterVariant)
        "Codons": "", # FILL? Should be like aGt/act (NucleotidsReference/NucleotidsVariant)
        "Existing_variation": rs.replace(";", "&"),
        "DISTANCE": ann_dict["Distance"],
        "STRAND": "-1",
        "FLAGS": ann_dict["Comments"],
        "VARIANT_CLASS": variant_class,
        "SYMBOL_SOURCE": "HGNC",
        "HGNC_ID": "", # FILL? Should be HGNC.XXX where XXX is the unique hugo id, no symbol
        "CANONICAL": "", # FILL? Should be an unknown flag indicating if the transcript is denoted as the canonical transcript for this gene
        "MANE": "", # ???
        "TSL": "", # FILL? Transcript support level
        "APPRIS": "", # FILL? APPRIS score
        "CCDS": "",
        "ENSP": "",
        "SWISSPROT": "",
        "TREMBL": "",
        "UNIPARC": "",
        "GENE_PHENO": "",
        "SIFT": "",
        "PolyPhen": "",
        "DOMAINS": "",
        "miRNA": "",
        "HGVS_OFFSET": "",
        "AF": "",
        "AFR_AF": "",
        "AMR_AF": "",
        "EAS_AF": "",
        "EUR_AF": "",
        "SAS_AF": "",
        "AA_AF": "",
        "EA_AF": "",
        "gnomAD_AF": "",
        "gnomAD_AFR_AF": "",
        "gnomAD_AMR_AF": "",
        "gnomAD_ASJ_AF": "",
        "gnomAD_EAS_AF": "",
        "gnomAD_FIN_AF": "",
        "gnomAD_NFE_AF": "",
        "gnomAD_OTH_AF": "",
        "gnomAD_SAS_AF": "",
        "MAX_AF": "",
        "MAX_AF_POPS": "",
        "CLIN_SIG": "",
        "SOMATIC": "",
        "PHENO": "",
        "PUBMED": "",
        "MOTIF_NAME": "",
        "MOTIF_POS": "",
        "HIGH_INF_POS": "",
        "MOTIF_SCORE_CHANGE": "",        
        "LoFtool": "",
    }

with (
    open(file=snakemake.input["vcf"], mode="r") as vcfinstream,
    open(file=snakemake.output["vcf"], mode="w") as vcfoutstream
):
    codons = {
        ""
    }
    for line in vcfinstream:
        if line.startswith("##"):
            pass
        elif line.startswith("#"):
            vcfoutstream.write(headers)
        else:
            chomp = line.split("\t")
            info = chomp[8]
            vep = snpeff2vep(get_snpeff_ann(chomp[8]))
            if vep is not None:
                info += ";" + vep
                chomp[8] = info
            line = "\t".join(chomp)
        

        vcfoutstream.write(line)

