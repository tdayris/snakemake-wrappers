rule pandas_filter_tsv:
   input:
      table="snpeff_snpsift/snpsift/extractFields/{sample}.tsv"
   output:
      table="pandas/filter/{how}/dp{dp}/{sample}.dp{dp}.tsv",
      xlsx="pandas/filter/{how}/dp{dp}/{sample}.dp{dp}.xlsx"
   threads: 1
   resources:
      mem_mb=lambda wildcards, attempt: attempt * 1024 * 8,
      time_min=lambda wildcards, attempt: attempt * 15,
      tmpdir="tmp"
   log:
      "logs/pandas_filter/{sample}/{how}.{dp}.log"
   params:
      new_cols = lambda wildcards: [
         ["Mutect2_Allele_Frequency", "=", f"{wildcards.sample}_tumor_AF"]
      ],
      prefixes = [
         ["Chromosome", "chr"]
      ],
      keep_column = lambda wildcards: [
            "Chromosome",
            "Start_Position",
            "Variant_ID",
            "SYMBOL",
            "Reference_Allele",
            "Tumor_Seq_Allele1",
            "Mutect2_Read_depth",
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
            "dbSNP_Pubmed",
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
            "Center",
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
            #[f"{wildcards.sample}_normal_AD_allele2", "<=", float(wildcards.dp)],
            #[f"{wildcards.sample}_tumor_AD_allele2", ">=", float(wildcards.dp)],
            [f"{wildcards.sample}_tumor_DP", ">=", float(wildcards.dp)],
            #["dbNSFP_ExAC_AlleleFrequency", "<=", 0.05],
            ["VarOcc", "<=", len(design["Sample_id"].tolist()) - 1]
      ] if config.get("ANMO", False) is True else [],
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
            #["Filter", "VarOccAbove105"],
            #["Filter", f"DepthBelow{wildcards.filter}X"],
            #["Filter", "AFBelow40pct"]
      ],
      contains = lambda wildcards: [
            [
                "Filter", 
                (
                    "ExistsInCancerGeneCensus" if wildcards.how == "census_only" else (
                        "ExistsInOncoKB" if wildcards.how == "oncokb_only" else ""
                    )
                )
            ]
      ],
      drop_duplicated_lines=True
   wrapper:
        "bio/pandas/filter_table"