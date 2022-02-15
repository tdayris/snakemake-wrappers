from pathlib import Path

tsv_dir_path = Path("snpsift/extractFields")
tsv_paths = [p.name for p in tsv_dir_path.iterdir() if p.name.endswith("tsv")]
sample_names = [tsv[:-len(".tsv")] for tsv in tsv_paths]

how_list = ["census_only", "oncokb_only", "all"]
dp_list = list(map(str, [5, 10, 40, 60]))

non_synonymous = [
    "Missense_Mutation",   "Nonsense_Mutation",    "In_Frame_Del",
    "3_prime_UTR_variant", "5_prime_UTR_variant",   "Frame_Shift_Del",
    "Splice_Site",         "Frame_Shift_Ins",       "In_Frame_Ins"
]
color_named_vector = [
    "#006400", "#dc143c", "#87cefa", "#1e90ff", "#0014a8",
    "#800080", "#483d8b", "#ca2c92", "#8e3a59"
]
gene_list = [
    "TP53",  "CDKN2A", "EGFR",   "FAT1", "KRAS",
    "MUC16", "NOCH1",  "PIK3CA", "LRP1B"
]

wildcard_constraints:
    how=r"|".join(how_list),
    dp=r"|".join(dp_list),
    sample=r"|".join(sample_names)


rule target:
   input:
      expand("results_to_upload/{how}/dp{dp}/{sample}.dp{dp}.tsv", how=how_list, sample=sample_names, dp=dp_list)


rule pandas_filter_tsv:
   input:
      table="snpsift/extractFields/{sample}.tsv"
   output:
      table="results_to_upload/{how}/dp{dp}/{sample}.dp{dp}.tsv"
   threads: 1
   resources:
      mem_mb=lambda wildcards, attempt: attempt * 1024 * 8,
      time_min=lambda wildcards, attempt: attempt * 15,
      tmpdir="tmp"
   log:
      "logs/pandas_filter/{sample}/{how}.{dp}.log"
   params:
      new_cols = lambda wildcards: [
         ["Mutect2_Allele_Frequency", "=", f"{wildcards.sample}_tumor_AF"],
         ["Chromosome", "=", "CHROM"],
         ["Reference_Allele", "=", "REF"],
         ["Variant_Type", "=", "VARTYPE"],
         ["Variant_ID", "=", "ID"],
         ["Filter", "=", "FILTER"],
         ["Variant_Classification", "=", "ANN[*].EFFECT"],
         ["Transcript_ID", "=", "ANN[*].FEATUREID"],
         ["Gene", "=", "ANN[*].GENEID"],
         ["dbNSFP_ExAC_AlleleFrequency", "=", "dbNSFP_ExAC_AF"],
         ["SYMBOL", "=", "ANN[*].GENE"],
         ["Mutect2_Read_depth", "=", f"{wildcards.sample}_tumor_DP"],
         ["Kaviar_Allele_Frequency", "=", "Kaviar_AF"],
         ["Hugo_Symbol", "=", "ANN[*].GENE"],
         ["HGVSc", "=", "ANN[*].HGVS_C"],
         ["Feature", "=", "ANN[*].FEATURE"],
         ["dbSNP_Clinical_Diagnostic_Assay", "=", "dbSNP_CDA"],
         ["MSigDb_Pathways", "=", "MSigDb"],
         ["BIOTYPE", "=", "ANN[*].BIOTYPE"],
         ["IMPACT", "=", "ANN[*].IMPACT"],
         ["HGVSp", "=", "ANN[*].HGVS_P"]
      ],
      prefixes = [
         ["Chromosome", "chr"]
      ],
      keep_column = lambda wildcards: [
            "Chromosome",
            "POS",
            "Variant_ID",
            "SYMBOL",
            "Reference_Allele",
            "Tumor_Seq_Allele1",
            "Tumor_Seq_Allele2",
            "HGVSc",
            "HGVSp",
            "VarOcc",
            "BIOTYPE",
            "IMPACT",
            "Variant_Classification",
            "Mutect2_Read_depth",
            "Mutect2_Allele_Frequency",
            "Variant_Type",
            "Tumor_Sample_Barcode",
            "dbSNP_Clinical_Diagnostic_Assay",
            "dbNSFP_ExAC_AlleleFrequency",
            "Kaviar_Allele_Frequency",
            "MSigDb_Pathways",
            "dbNSFP_Polyphen2_HDIV_pred",
            "dbNSFP_SIFT_pred",
            "dbNSFP_FATHMM_pred",
            "dbNSFP_ClinPred_score",
            "dbNSFP_MutationTaster_pred",
            'dbNSFP_MutationAssessor_pred',
            "dbNSFP_MetaLR_pred",
            "dbNSFP_LIST_S2_pred",
            "dbNSFP_LRT_pred",
            "dbNSFP_BayesDel_noAF_pred",
            "dbNSFP_Aloft_pred",
            "Transcript_ID",
            "Gene",
            "Feature",
            "Filter",
            "Hugo_Symbol",
            "End_Position",
            "CancerGeneCensus_Gene_Symbol",
            "CancerGeneCensus_Name",
            "CancerGeneCensus_Entrez_GeneId",
            "CancerGeneCensus_Tier",
            "CancerGeneCensus_Somatic",
            "CancerGeneCensus_Germline",
            "CancerGeneCensus_Tumour_TypesSomatic",
            "OncoKB_Hugo_Symbol",
            "OncoKB_Entrez_Gene_ID",
            "OncoKB_GRCh37_Isoform",
            "OncoKB_GRCh37_RefSeq",
            "OncoKB_GRCh38_RefSeq",
            "OncoKB_OncoKB_Annotated",
            "OncoKB_Is_Oncogene",
            "OncoKB_Is_Tumor_Suppressor_Gene",
            "OncoKB_MSK_IMPACT",
            "OncoKB_MSK_HEME",
            f"{wildcards.sample}_normal_Reference_Allele",
            f"{wildcards.sample}_normal_Seq_Allele1",
            f"{wildcards.sample}_normal_Seq_Allele2",
            f"{wildcards.sample}_normal_DP",
            f"{wildcards.sample}_normal_AD_allele2",
            f"{wildcards.sample}_normal_AF",
            f"{wildcards.sample}_tumor_Reference_Allele",
            f"{wildcards.sample}_tumor_Seq_Allele1",
            f"{wildcards.sample}_tumor_Seq_Allele2",
            f"{wildcards.sample}_tumor_DP",
            f"{wildcards.sample}_tumor_AD_allele1",
            f"{wildcards.sample}_tumor_AD_allele2",
            f"{wildcards.sample}_tumor_AF",
      ],
      convert_cols_type = lambda wildcards: {
            "Mutect2_Allele_Frequency": "float",
            "Mutect2_Read_depth": "int",
            f"{wildcards.sample}_tumor_AF": "float",
            f"{wildcards.sample}_normal_AF": "float",
            #f"{wildcards.sample}_normal_AD_allele2": "float",
            f"{wildcards.sample}_tumor_DP": "int",
            #f"{wildcards.sample}_tumor_AD_allele2": "float",
      },
      filters = lambda wildcards: [
            # ["Variant_Classification", "!=", "downstream_gene_variant"],
            # ["Variant_Classification", "!=", "intergenic_region"],
            # ["Variant_Classification", "!=", "synonymous_variant"],
            # ["Variant_Classification", "!=", "non_coding_transcript_exon_variant"],
            # ["Variant_Classification", "!=", "upstream_gene_variant"],
            # ["Variant_Classification", "!=", "splice_region_variant&synonymous_variant"],
            # ["Variant_Classification", "!=", "non_coding_transcript_variant"],
            # ["Variant_Classification", "!=", "intron_variant"],
            ["Variant_Classification", "!=", "Synonymous_Variant"],
            #["Mutect2_Read_depth", ">=", int(wildcards.filter)],
            #['Mutect2_Allele_Frequency', ">=", 0.1],
            [f"{wildcards.sample}_tumor_AF", ">=", 0.1],
            [f"{wildcards.sample}_normal_AF", "<=", 0.1],
            ##[f"{wildcards.sample}_normal_AD_allele2", "<=", float(wildcards.dp)],
            ##[f"{wildcards.sample}_tumor_AD_allele2", ">=", float(wildcards.dp)],
            [f"{wildcards.sample}_tumor_DP", ">=", float(wildcards.dp)],
            #["dbNSFP_ExAC_AlleleFrequency", "<=", 0.05],
            ["VarOcc", "<=", float(len(sample_names) - 1)]
      ],
      not_contains = lambda wildcards: [
            ["Filter", "germline"],
            ["Filter", "multiallelic"],
            #"Filter": "fragment",
            #"Filter": "contamination",
            #"Filter": "weak_evidence",
            #"Filter": "slippage",
            #"Filter": "strand_bias",
            #"Filter": "map_qual",
            ["Filter", "haplotype"],
            #"Filter": "base_qual",
            #["Filter", "DPBelow5"],
            ["Filter", "AboveFisherStrandBias"],
            ["Filter", "AboveStrandOddsRatio"],
            ["Filter", "BelowMQRankSum"],
            ["Filter", "BelowReadPosRankSum"],
            ##["Filter", "VarOccAbove105"],
            #["Filter", f"DepthBelow{wildcards.filter}X"],
            #["Filter", "AFBelow40pct"]
      ],
      contains = lambda wildcards: [
            ["Filter", (
               "ExistsInCanceGeneCensus" if wildcards.how == "census_only" else (
                  "ExistsInOncoKB" if wildcards.how == "oncokb_only" else ""))]
      ],
      drop_duplicated_lines=True
   wrapper:
        "bio/pandas/filter_table"
