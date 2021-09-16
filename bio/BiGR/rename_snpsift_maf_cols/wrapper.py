#!/usr/bin/python3.8
# -*- coding: utf-8 -*-

"""
Rename splitted vcf columns
"""

import pandas
import logging

logging.basicConfig(
    filename=snakemake.log[0],
    filemode="w",
    level=logging.DEBUG
)

def get_variant_type(ref: str, alt: str) -> str:
    """
    Return variant type given the reference and the alternative alleles
    """
    if ref == ".": # Nothing in REF, added in ALT
        return "INS"
    if alt == ".":
        return "DEL" # Something in REF, nothing in ALT

    if len(ref) == 1:
        if len(alt) == 1:
            return "SNP" # REF and ALT only have one nucleotide
        return "ONP" # REF had less nucleotides than ALT

    if len(ref) == 2:
        if len(alt) == 2:
            return "DNP" # REF and ALT have two nucleotides
        return "ONP" # REF has two nucleotides, ALT have unknown

    if len(ref) == 3:
        if len(alt) == 3:
            return "TNP" # REF and ALT have three nucleotides
        return "ONP" # REF has three nucleotides, ALT have unknown
    return "ONP" # More than 3 nucleotides in ref, unknown in ALT


headers = {
    "CHROM": "Chromosome",
    "POS": "Start_Position",
    "ID": "Variant_ID",
    "REF3": "Reference_Allele",
    "ALT": "Tumor_Seq_Allele1",
    "FILTER": "Filter",
    "FORMAT_AD": "Allelic_depths",
    "FORMAT_AF": "Allele_fractions",
    "FORMAT_DP": "Pileup_Read_Depth",
    "FORMAT_F1R2": "F1R2_orientation",
    "FORMAT_F2R1": "F2R1_orientation",
    "FORMAT_GQ": "Genotype_quality",
    "FORMAT_GT": "Genotype",
    "FORMAT_PGT": "Physical_phasing_haplotype",
    "FORMAT_PID": "Physical_phasing_ID",
    "FORMAT_PL": "Normalized_Phred_scaled_likelihoods",
    "FORMAT_PS": "Phasing_set",
    "FORMAT_SB": "Strand_bias_Fisher_component",
    "AS_SB_TABLE": "Allele_specific_Strand_bias",
    "AS_UNIQ_ALT_READ_COUNT": "Unique_alt_variant_count",
    "CONTQ": "Phred_scaled_qualities_non_contamination",
    "DP": "Mutect2_Read_depth",
    "ECNT": "Events_in_haplotype",
    "GERMQ": "Phred_scaled_quality_non_germline",
    "MBQ": "Median_base_quality",
    "MFRL": "Median_fragment_length",
    "MMQ": "Median_mapping_quality",
    "MPOS": "Median_distance_to_end_read",
    "NALOD": "Negative_log_10_odds_of_artifact",
    "NCount": "Pileup_N_base_depth",
    "NLOD": "Normal_log_10_likelihood_ploïdy",
    "OCM": "Non_matching_original_alignment",
    "PON": "Panel_of_Normal",
    "POPAF": "negative_log_10_population_allele_frequencies",
    "ROQ": "Phred_scaled_qualities_not_orientation_bias",
    "RPA": "Number_tandem_repetition",
    "RU": "Repeat_unit",
    "SEQQ": "Phred_scaled_quality_not_sequencing_error",
    "STR": "Is_short_tandem_repeat",
    "STRANDQ": "Phred_scaled_quality_not_strand_bias",
    "STRQ": "Phred_scaled_quality_not_tandem_polymerase_error",
    "TLOD": "Log_10_likelihood_variant_exists",
    "ANN[*].ALLELE": "SnpEff_Genotype",
    "ANN[*].EFFECT": "Variant_Classification",
    "ANN[*].IMPACT": "IMPACT",
    "ANN[*].GENE": "Hugo_Symbol",
    "ANN[*].GENEID": "Gene",
    "ANN[*].FEATURE": "Feature",
    "ANN[*].FEATUREID": "Transcript_ID",
    "ANN[*].BIOTYPE": "BIOTYPE",
    "ANN[*].RANK": "Exon_Intron_rank",
    "ANN[*].HGVS_C": "HGVSc",
    "ANN[*].HGVS_P": "HGVSp",
    "ANN[*].CDNA_POS": "cDNA_position",
    "ANN[*].CDNA_LENANN[*].CDS_POS": "CDS_position",
    "ANN[*].CDS_LEN": "CDS_length",
    "ANN[*].AA_POS": "Protein_position",
    "ANN[*].AA_LEN": "Protein_length",
    "ANN[*].DISTANCE": "DISTANCE",
    "ANN[*].ERRORS": "Errors",
    "LOF": "Loss_of_Function_SnpEff",
    "NMD": "Nonsense_mediated_decay",
    "VARTYPE": "Variant_types",
    "SNP": "Is_SNP",
    "MNP": "Is_MNP",
    "INS": "Is_INS",
    "DEL": "Is_DEL",
    "MIXED": "Is_MIXED_Polymorphisms",
    "HOM": "Is_Homozygous",
    "HET": "Is_Heterozygous",
    "MSigDb": "MSigDb_Pathway",
    "AC": "Allele_Count",
    "AF": "Allele_Frequency",
    "AN": "Allele_Number",
    "DS": "Database_Source",
    "END": "End_Position",
    "CDA": "Clinical_Diagnostic_Assay",
    "OTH": "Orthologous_Variants_in_NCBI",
    "S3D": "Has_known_3D_structure",
    "WTD": "Is_Withdrawn_by_submitter",
    "dbSNPBuildID": "dbSNP_id",
    "SLO": "SubmitterLinkOut",
    "NSF": "Has_Non_Synonymous_Frameshift",
    "R3": "3p_variant",
    "R5": "5p_variant",
    "NSN": "Has_Nonsense_Changes_Stop",
    "NSM": "Has_Nonsense_Changes_Peptide",
    "G5A": "Less_5pct_population",
    "COMMON": "Is_common_variant",
    "RS": "RS_id_dbSNP",
    "RV": "Variant_is_reversed",
    "TPA": "Third_party_annotation",
    "CFL": "Has_assembly_conflict",
    "GNO": "Has_Available_genotype",
    "VLD": "Is_Validated",
    "ASP": "Is_Assembly_Specific",
    "ASS": "Acceptor_splice_variant",
    "REF": "Variant_identical_to_reference",
    "U3": "3p_UTR_variant",
    "U5": "5p_UTR_variant",
    "TOPMED": "AF_in_TopMed_db",
    "WGT": "Weight",  # Weight, 00 - unmapped, 1 - weight 1, 2 - weight 2, 3 - weight 3 or more
    "MTP": "Microattribution_TPA",
    "LSD": "Locus_specific_database",
    "NOC": "Contig_not_in_variant_list",
    "DSS": "Donor_splice_varianP",
    "SYN": "Synonymous_variant",
    "KGPhase3": "1000genome_p3",
    "CAF": "Alleles_in_1000g",
    "VC": "Variant_class",
    "MUT": "Is_cited",
    "KGPhase1": "1000genome_p1",
    "VP": "Variant_property",
    "SAO": "Variant_Origin", # Origin: 0 - unspecified, 1 - Germline, 2 - Somatic, 3 - Both"
    "GENEINFO": "Gene_info",
    "INT": "Is_Intronic",
    "G5": "MAF_Over_5pct",
    "OM": "Has_OMIM_OMIA",
    "PMC": "Pubmed",
    "SSR": "Suspect_reason_code", # 0 - unspecified, 1 - Paralog, 2 - byEST, 4 - oldAlign, 8 - Para_EST, 16 - 1kg_failed, 1024 - other"
    "RSPOS": "Position_dbSNP",
    "HD": "High_density_genotyping_kit_dbGaP",
    "PM": "Variant_Precious",
    "QD": "Variant_call_confidence",
    "MQRankSum": "MQRankSum",
    "RAW_MQ": "Raw_Mapping_Quality_Database",
    "FS": "FisherStrand_bias_dbSNP",
    "POSITIVE_TRAIN_SITE": "Used_by_GATK_VQSR_training_pos",
    "culprit": "Worst_annotation_VQSR",
    "NEGATIVE_TRAIN_SITE": "Used_by_GATK_VQSR_training_neg",
    "variant_type": "dbSNP_Variant_Type",
    "ReadPosRankSum": "Position_bias_dbSNP",
    "SOR": "SymetricOddRatio_strand_bias_dnSNP",
    "InbreedingCoeff": "Inbreeding_Coefficient_Hardy_Weinberg",
    "vep": "Vep_Annotation",
    "MQ": "Mapping_Quality_dnSNP",
    "n_alt_alleles": "Alternative_Alleles_dbSNP",
    "lcr": "Low_Complexity_Region",
    "AS_VQSLOD": "GATK_VQSR_log_Odds_Ratio",
    "faf95_nfe": "Filtered_Allele_Frequency_95pct_Non_Finnish_Europeans",
    "AC_asj": "Alternative_Count_Ashkenazi_Jewish_ancestry",
    "nhomalt_ami_female": "Homozygous_individuals_Amish_Ancestry_female",
    "AN_fin_male": "Allele_Number_Finnish_ancestry_male",
    "nhomalt_eas": "Homozygous_individuals_East_Asian_ancestry",
    "AN_nfe_male": "Allele_Number_Non_Finnish_Europeans_male",
    "AC_male": "Alternative_Count_male",
    "AF_oth": "Allele_Frequency_Other_Ancestry",
    "nhomalt_ami_male": "Homozygous_individuals_Amish_Ancestry_male",
    "AC_nfe_male": "Alternative_Count_Non_Finnish_Europeans_male",
    "AN_sas_male": "Allele_Number_South_Asian_ancestry_male",
    "AF_nfe_male": "Allele_Frequency_Non_Finnish_Europeans_male",
    "AF_afr": "Allele_Frequency_African_Ancestry",
    "AC_amr_male": "Alternative_Count_Latino_Ancestry_male",
    "nhomalt_amr_female": "Homozygous_individuals_Latino_Ancestry_female",
    "nhomalt_oth": "Homozygous_individuals_Other_Ancestry",
    "nhomalt_oth_male": "Homozygous_individuals_Other_Ancestry_male",
    "nhomalt_eas_female": "Homozygous_individuals_East_Asian_ancestry_female",
    "AN_oth_male": "Allele_Number_Other_Ancestry_male",
    "AC_nfe_female": "Alternative_Count_Non_Finnish_Europeans_female",
    "nhomalt_nfe_female": "Homozygous_individuals_Non_Finnish_Europeans_female",
    "faf95_sas": "Filtered_Allele_Frequency_95pct_South_Asian_ancestry",
    "nhomalt_nfe_male": "Homozygous_individuals_Non_Finnish_Europeans_male",
    "AN_ami_female": "Allele_Number_Amish_Ancestry_female",
    "AC_oth": "Alternative_Count_Other_Ancestry",
    "AC_female": "Alternative_Count_female",
    "faf95_afr": "Filtered_Allele_Frequency_95pct_African_Ancestry",
    "nhomalt_fin": "Homozygous_individuals_Finnish_ancestry",
    "AC_ami": "Alternative_Count_Amish_Ancestry",
    "AC_eas": "Alternative_Count_East_Asian_ancestry",
    "AF_asj": "Allele_Frequency_Ashkenazi_Jewish_ancestry",
    "AF_amr": "Allele_Frequency_Latino_Ancestry",
    "AF_oth_male": "Allele_Frequency_Other_Ancestry_male",
    "AF_fin_female": "Allele_Frequency_Finnish_ancestry_female",
    "AC_afr_male": "Alternative_Count_African_Ancestry_male",
    "AN_female": "Allele_Number_female",
    "AC_eas_female": "Alternative_Count_East_Asian_ancestry_female",
    "AN_ami": "Allele_Number_Amish_Ancestry",
    "AN_afr_male": "Allele_Number_African_Ancestry_male",
    "faf99_sas": "Filtered_Allele_Frequency_99pct_South_Asian_ancestry",
    "nhomalt_sas_male": "Homozygous_individuals_South_Asian_ancestry_male",
    "AF_male": "Allele_Frequency_male",
    "nhomalt_afr": "Homozygous_individuals_African_Ancestry",
    "AC_asj_male": "Alternative_Count_Ashkenazi_Jewish_ancestry_male",
    "AC_raw": "Alternative_Count_raw",
    "AN_oth": "Allele_Number_Other_Ancestry",
    "AN_oth_female": "Allele_Number_Other_Ancestry_female",
    "AN_amr_female": "Allele_Number_Latino_Ancestry_female",
    "AN_amr_male": "Allele_Number_Latino_Ancestry_male",
    "AF_amr_female": "Allele_Frequency_Latino_Ancestry_female",
    "AF_sas": "Allele_Frequency_South_Asian_ancestry",
    "AF_eas_female": "Allele_Frequency_East_Asian_ancestry_female",
    "AN_asj_male": "Allele_Number_Ashkenazi_Jewish_ancestry_male",
    "faf99_amr": "Filtered_Allele_Frequency_99pct_Latino_Ancestry",
    "nhomalt_oth_female": "Homozygous_individuals_Other_Ancestry_female",
    "AN_eas": "Allele_Number_East_Asian_ancestry",
    "AC_sas_female": "Alternative_Count_South_Asian_ancestry_female",
    "AF_nfe_female": "Allele_Frequency_Non_Finnish_Europeans_female",
    "AN_asj_female": "Allele_Number_Ashkenazi_Jewish_ancestry_female",
    "nhomalt_male": "Homozygous_individuals_male",
    "AN_asj": "Allele_Number_Ashkenazi_Jewish_ancestry",
    "AN_nfe_female": "Allele_Number_Non_Finnish_Europeans_female",
    "AN_ami_male": "Allele_Number_Amish_Ancestry_male",
    "AC_sas_male": "Alternative_Count_South_Asian_ancestry_male",
    "nhomalt_afr_female": "Homozygous_individuals_African_Ancestry_female",
    "AC_nfe": "Alternative_Count_Non_Finnish_Europeans",
    "AF_ami_male": "Allele_Frequency_Amish_Ancestry_male",
    "AN_amr": "Allele_Number_Latino_Ancestry",
    "AF_nfe": "Allele_Frequency_Non_Finnish_Europeans",
    "AN_fin_female": "Allele_Number_Finnish_ancestry_female",
    "AN_eas_female": "Allele_Number_East_Asian_ancestry_female",
    "AC_ami_female": "Alternative_Count_Amish_Ancestry_female",
    "AN_eas_male": "Allele_Number_East_Asian_ancestry_male",
    "AF_eas": "Allele_Frequency_East_Asian_ancestry",
    "faf99_eas": "Filtered_Allele_Frequency_99pct_East_Asian_ancestry",
    "nhomalt_fin_male": "Homozygous_individuals_Finnish_ancestry_male",
    "AN_sas": "Allele_Number_South_Asian_ancestry",
    "nhomalt_raw": "Homozygous_individuals_raw",
    "nhomalt_sas_female": "Homozygous_individuals_South_Asian_ancestry_female",
    "nhomalt_asj": "Homozygous_individuals_Ashkenazi_Jewish_ancestry",
    "AF_ami": "Allele_Frequency_Amish_Ancestry",
    "faf95_eas": "Filtered_Allele_Frequency_95pct_East_Asian_ancestry",
    "nhomalt_nfe": "Homozygous_individuals_Non_Finnish_Europeans",
    "nhomalt_afr_male": "Homozygous_individuals_African_Ancestry_male",
    "AC_fin": "Alternative_Count_Finnish_ancestry",
    "AF_oth_female": "Allele_Frequency_Other_Ancestry_female",
    "AF_asj_male": "Allele_Frequency_Ashkenazi_Jewish_ancestry_male",
    "AC_eas_male": "Alternative_Count_East_Asian_ancestry_male",
    "AC_amr": "Alternative_Count_Latino_Ancestry",
    "AF_sas_male": "Allele_Frequency_South_Asian_ancestry_male",
    "AC_oth_male": "Alternative_Count_Other_Ancestry_male",
    "AN_male": "Allele_Number_male",
    "faf95_adj": "Filtered_Allele_Frequency_95pct_adj",
    "nhomalt_eas_male": "Homozygous_individuals_East_Asian_ancestry_male",
    "AF_fin_male": "Allele_Frequency_Finnish_ancestry_male",
    "AC_sas": "Alternative_Count_South_Asian_ancestry",
    "AC_afr": "Alternative_Count_African_Ancestry",
    "AF_ami_female": "Allele_Frequency_Amish_Ancestry_female",
    "faf95_amr": "Filtered_Allele_Frequency_95pct_Latino_Ancestry",
    "AF_afr_female": "Allele_Frequency_African_Ancestry_female",
    "AF_fin": "Allele_Frequency_Finnish_ancestry",
    "AF_eas_male": "Allele_Frequency_East_Asian_ancestry_male",
    "faf99_nfe": "Filtered_Allele_Frequency_99pct_Non_Finnish_Europeans",
    "AC_oth_female": "Alternative_Count_Other_Ancestry_female",
    "AF_female": "Allele_Frequency_female",
    "AF_amr_male": "Allele_Frequency_Latino_Ancestry_male",
    "nhomalt_fin_female": "Homozygous_individuals_Finnish_ancestry_female",
    "AN_raw": "Allele_Number_raw",
    "nhomalt_amr": "Homozygous_individuals_Latino_Ancestry",
    "nhomalt_amr_male": "Homozygous_individuals_Latino_Ancestry_male",
    "AC_ami_male": "Alternative_Count_Amish_Ancestry_male",
    "nhomalt_asj_female": "Homozygous_individuals_Ashkenazi_Jewish_ancestry_female",
    "faf99_adj": "Filtered_Allele_Frequency_99pct_adj",
    "AN_afr_female": "Allele_Number_African_Ancestry_female",
    "AC_amr_female": "Alternative_Count_Latino_Ancestry_female",
    "AF_asj_female": "Allele_Frequency_Ashkenazi_Jewish_ancestry_female",
    "AN_sas_female": "Allele_Number_South_Asian_ancestry_female",
    "AF_sas_female": "Allele_Frequency_South_Asian_ancestry_female",
    "nhomalt_ami": "Homozygous_individuals_Amish_Ancestry",
    "nhomalt_sas": "Homozygous_individuals_South_Asian_ancestry",
    "nhomalt": "Homozygous_individuals",
    "AC_fin_male": "Alternative_Count_Finnish_ancestry_male",
    "nhomalt_asj_male": "Homozygous_individuals_Ashkenazi_Jewish_ancestry_male",
    "AC_asj_female": "Alternative_Count_Ashkenazi_Jewish_ancestry_female",
    "AF_afr_male": "Allele_Frequency_African_Ancestry_male",
    "AN_fin": "Allele_Number_Finnish_ancestry",
    "AC_fin_female": "Alternative_Count_Finnish_ancestry_female",
    "nhomalt_female": "Homozygous_individuals_female",
    "AN_afr": "Allele_Number_African_Ancestry",
    "AC_afr_female": "Alternative_Count_African_Ancestry_female",
    "dbNSFP_ExAC_NFE_AF": "ExAC_AF_NFE",
    "dbNSFP_ExAC_SAS_AF": "ExAC_AF_SAS",
    "dbNSFP_GERP___RS": "GERPpp_RS",
    "dbNSFP_GERP___NR": "GERPpp_NR",
    "dbNSFP_ExAC_Adj_AC": "ExAC_Adjusted_Allele_Count",
    "dbNSFP_ExAC_Adj_AF": "ExAC_Adjusted_Allele_Frequency",
    "dbNSFP_ExAC_SAS_AC": "ExAC_AC_SAS",
    "dbNSFP_1000Gp3_AMR_AF": "AMR_MAF",
    "dbNSFP_1000Gp3_AMR_AC": "AMR_MAF_AC",
    "dbNSFP_MetaSVM_pred": "MetaSNV_prediction_score",
    "dbNSFP_1000Gp3_EAS_AC": "EAS_MAF_AC",
    "dbNSFP_Interpro_domain": "Interpro_domain",
    "dbNSFP_FATHMM_pred": "FATHMM_prediction_score",
    "dbNSFP_ExAC_AFR_AF": "ExAC_AF_AFR",
    "dbNSFP_ExAC_AFR_AC": "ExAC_AC_AFR",
    "dbNSFP_1000Gp3_AC": "1000g_general_allele_count",
    "dbNSFP_1000Gp3_AF": "GMAF",
    "dbNSFP_1000Gp3_EAS_AF": "EAS_MAF",
    "dbNSFP_ExAC_AF": 'ExAC_AF',
    "dbNSFP_Uniprot_acc": "TREMBL",
    "dbNSFP_ExAC_AC": "ExAC_AC",
    "dbNSFP_LRT_pred": "LRT_prediction_score",
    "dbNSFP_PROVEAN_pred": "Provean_prediction_score",
    "dbNSFP_ExAC_FIN_AC": "ExAC_AC_FIN",
    "dbNSFP_phastCons100way_vertebrate": "phastCons100way_vertebrate_prediction_score",
    "dbNSFP_ExAC_FIN_AF": "ExAC_AF_FIN",
    "dbNSFP_CADD_phred": "CADD_phred_score",
    "dbNSFP_1000Gp3_EUR_AC": "EUR_MAF_AC",
    "dbNSFP_Polyphen2_HDIV_pred": "PolyPhen",
    "dbNSFP_ESP6500_EA_AC": "Allele_count_Europe_American_NHLBIGO",
    "dbNSFP_1000Gp3_AFR_AC": 'AFR_MAF_AC',
    "dbNSFP_1000Gp3_EUR_AF": "EUR_MAF",
    "dbNSFP_ExAC_AMR_AF": "ExAC_AF_AMR",
    "dbNSFP_1000Gp3_AFR_AF": 'AFR_MAF',
    'dbNSFP_MutationTaster_pred': "Mutation_taster_prediction_score",
    "dbNSFP_MutationAssessor_pred": "Mutation_assessor_prediction_score",
    "dbNSFP_ESP6500_AA_AF": "Allele_frequency_African_American_NHLBIGO",
    "dbNSFP_Polyphen2_HVAR_pred": "Polyphen_HVAR_prediction_score",
    "dbNSFP_ExAC_AMR_AC": "ExAC_AC_AMR",
    "dbNSFP_ExAC_NFE_AC": "ExAC_AC_NFE",
    "dbNSFP_SIFT_pred": "SIFT",
    "dbNSFP_1000Gp3_SAS_AC": "SAS_MAF_AC",
    "dbNSFP_ExAC_EAS_AC": "ExAC_AC_EAS",
    "dbNSFP_1000Gp3_SAS_AF": "SAS_MAF",
    "dbNSFP_ExAC_EAS_AF": "ExAC_AF_EAS",
    "dbNSFP_ESP6500_EA_AF": "Allele_frequency_Europe_American_NHLBIGO",
    "dbNSFP_ESP6500_AA_AC": "Allele_count_African_American_NHLBIGO",
}

# Gather parameters
add_cols = snakemake.params.get("add_cols", True)
ncbi_build = snakemake.params.get("ncbi_build", "GRCh38")
center = snakemake.params.get("center", "GustaveRoussy")
#caller = snakemake.params.get("caller", "mutect2")
logging.info("Parameters retrieved")

# Load user's data
variants = pandas.read_csv(
    snakemake.input["tsv"],
    sep="\t",
    header=0,
    index_col=None
)
logging.info("Variants loaded in memory")

# Replace header names
new_header = []
# translation_table = (
#     headers_mutect2
#     if caller.lower() == "mutect2"
#     else None
# )
translation_table = headers

for idx, colname in enumerate(variants.columns.tolist()):
    new_colname = translation_table.get(colname, colname)
    if idx == 3 and colname == "REF":
        new_header.append("Reference_Allele")
    else:
        new_header.append(new_colname if new_colname is not None else colname)
variants.columns = new_header
logging.info("New header defined")
logging.debug(variants.columns.tolist())

# Add new columns on demand
if add_cols is True:
    if "Center" not in variants.columns:
        logging.debug("Adding Sequencing Center information")
        variants["Center"] = variants["Hugo_Symbol"]

    if "NCBI_Build" not in variants.columns:
        logging.debug("Adding NCBI build information")
        variants["NCBI_Build"] = variants["Hugo_Symbol"]

    if "End_Position" not in variants.columns:
        logging.debug("Adding End position")
        variants["End_Position"] = [
            start - (len(ref) - 1)
            for start, ref in zip(
                variants["Start_Position"], variants["Reference_Allele"]
            )
        ]

    if "Variant_Type" not in variants.columns:
        logging.debug("Adding variant type")
        logging.debug(variants.head())
        variants["Variant_Type"] = [
            get_variant_type(ref, alt)
            for ref, alt
            in zip(variants["Reference_Allele"], variants["Tumor_Seq_Allele1"])
        ]

    if "SYMBOL" not in variants.columns:
        logging.debug("Adding gene symbols")
        variants["SYMBOL"] = variants["Hugo_Symbol"]

    if "HGNC_ID" not in variants.columns:
        logging.debug("Adding hugo symbols")
        variants["HGNC_ID"] = variants["Hugo_Symbol"]

    if ("ExAC_AF_AFR" not in variants.columns) and ("Allele_Frequency_African_Ancestry" in variants.columns):
        logging.debug("Adding ExAC_AF_AFR symbols")
        variants["ExAC_AF_AFR"] = variants["Allele_Frequency_African_Ancestry"]

    if ("ExAC_AF_AMR" not in variants.columns) and ("Allele_Frequency_Latino_Ancestry" in variants.columns):
        logging.debug("Adding ExAC_AF_AMR symbols")
        variants["ExAC_AF_AMR"] = variants["Allele_Frequency_Latino_Ancestry"]

    if ("ExAC_AF_EAS" not in variants.columns) and ("Allele_Frequency_East_Asian_ancestry" in variants.columns):
        logging.debug("Adding ExAC_AF_EAS symbols")
        variants["ExAC_AF_EAS"] = variants["Allele_Frequency_East_Asian_ancestry"]

    if ("ExAC_AF_FIN" not in variants.columns) and ("Allele_Frequency_Finnish_ancestry" in variants.columns):
        logging.debug("Adding ExAC_AF_FIN symbols")
        variants["ExAC_AF_FIN"] = variants["Allele_Frequency_Finnish_ancestry"]

    if ("ExAC_AF_NFE" not in variants.columns) and ("Allele_Frequency_Non_Finnish_Europeans" in variants.columns):
        logging.debug("Adding ExAC_AF_NFE symbols")
        variants["ExAC_AF_NFE"] = variants["Allele_Frequency_Non_Finnish_Europeans"]

    if ("ExAC_AF_OTH" not in variants.columns) and ("Allele_Frequency_Other_Ancestry" in variants.columns):
        logging.debug("Adding ExAC_AF_OTH symbols")
        variants["ExAC_AF_OTH"] = variants["Allele_Frequency_Other_Ancestry"]

    if ("ExAC_AF_SAS" not in variants.columns) and ("Allele_Frequency_South_Asian_ancestry" in variants.columns):
        logging.debug("Adding ExAC_AF_SAS symbols")
        variants["ExAC_AF_SAS"] = variants["Allele_Frequency_South_Asian_ancestry"]

    if "vcf_region" not in variants.columns:
        logging.debug("Adding VCF regions")
        variants["vcf_region"] = [
            ":".join(map(str, [chr, pos, id, ref, alt]))
            for chr, pos, id, ref, alt in zip(
                variants["Chromosome"],
                variants["Start_Position"],
                variants["Variant_ID"],
                variants["Reference_Allele"],
                variants["Tumor_Seq_Allele1"]
            )
        ]

    if "Variant_Origin" in variants.columns:
        logging.debug("Translating Variant_Origin")
        translation = {
            "0": "unspecified",
            "0.0": "unspecified",
            "nan": "unspecified",
            "1": "germline",
            "1.0": "germline",
            "2": "somatic",
            "2.0": "somatic",
            "3": "both_germline_somatic",
            "3.0": "both_germline_somatic"
        }

        variants["Variant_Origin_Readable"] = [
            translation[str(i)] for i in variants["Variant_Origin"]
        ]

    if "Suspect_reason_code" in variants.columns:
        translation = {
            "0": "unspecified",
            "1": "Paralog",
            "2": "byEST",
            "4": "oldAlign",
            "8": "Para_EST",
            "16": "1kg_failed",
            "1024": "other",
            "0.0": "unspecified",
            "1.0": "Paralog",
            "2.0": "byEST",
            "4.0": "oldAlign",
            "8.0": "Para_EST",
            "16.0": "1kg_failed",
            "1024.0": "other",
            "nan": "unspecified"
        }

        variants["Suspect_reason_code_Readable"] = [
            translation[str(i)] for i in variants["Suspect_reason_code"]
        ]

    if "Tumor_Sample_Barcode" in snakemake.params.keys():
        variants["Tumor_Sample_Barcode"] = [
            snakemake.params["Tumor_Sample_Barcode"]
            for _ in variants["Start_Position"]
        ]

    if "Matched_Norm_Sample_Barcode" in snakemake.params.keys():
        variants["Matched_Norm_Sample_Barcode"] = [
            snakemake.params["Matched_Norm_Sample_Barcode"]
            for _ in variants["Start_Position"]
        ]

variants.to_csv(snakemake.output["tsv"], sep="\t", index=False)
