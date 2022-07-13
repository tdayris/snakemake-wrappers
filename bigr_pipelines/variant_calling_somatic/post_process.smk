from pathlib import Path

maf_path = Path("snpsift/fixed/")
vcf_list = [
    i.absolute()
    for i in maf_path.iterdir()
    if i.name.endswith("vcf.gz") and not i.name.startswith("complete")
    #if not i.name.startswith(("s087", "s083", "s250", "s157", "s198", "s094", "s184", "s249"))
]
sample_names = [i.name[:-len(".maf")] for i in vcf_list]

how_list = ["census_only", "all", "oncokb_only"]
dp_list = list(map(str, [5, 10, 40, 60]))
png_content = ["oncoplot", "summary", "titv"] #  , "signatures"]

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
    filter=r"|".join(dp_list),
    how=r"|".join(how_list),
    sample=r"|".join(sample_names),
    groups=r"alive|deceased",
    png=r"|".join(png_content),
    gene=r"|".join(gene_list)


ruleorder: gatk_variant_filtration > tabix_index
ruleorder: gatk_variant_filtration > pbgzip_compress


rule all:
    input:
        expand(
            "results_to_upload/{how}/dp{filter}/{sample}_dp{filter}.tsv",
            sample=sample_names,
            filter=dp_list,
            how=how_list
        ),
        expand(
            "results_to_upload/filtered_vcf/{sample}.vcf.gz{ext}",
            sample=sample_names,
            ext=["", ".tbi"]
        ),
        expand(
            "results_to_upload/{how}/{png}_dp{filter}.png",
            how=how_list,
            filter=dp_list,
            png=png_content
        ),
        expand(
            "results_to_upload/{how}/{gene}_dp{filter}.png",
            how=how_list,
            filter=dp_list,
            gene=gene_list
        )
    message:
        "Finishing ANMO Post process"


################
### MAFTools ###
################


rule maftools_lollipop:
    input:
        maf="tmp/{how}.dp{filter}.RDS"
    output:
        png="results_to_upload/{how}/{gene}_dp{filter}.png"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 5,
        time_min=lambda wildcards, attempt: attempt * 15,
        tmpdir="tmp"
    log:
        "logs/maftools/lollipop_plot/{how}_{filter}_{gene}.log"
    params:
        lollipop_plot_extra = "AACol = 'HGVSp', showMutationRate = TRUE"
    wrapper:
        "bio/maftools/lollipop_plot"


rule maftools_signatures:
    input:
        rds = "tmp/matrix_{how}.dp{filter}.RDS"
    output:
        png = "results_to_upload/{how}/signatures_dp{filter}.png",
        rds = temp("tmp/signatures/signatures_dp{filter}_{how}.RDS")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 5,
        time_min=lambda wildcards, attempt: attempt * 15,
        tmpdir="tmp"
    params:
        sig_extra = "n=2"
    wrapper:
        "bio/maftools/signatures"


rule maftools_trinucleotidematrix_hg38:
    input:
        rds = "tmp/{how}.dp{filter}.RDS"
    output:
        tsv = temp("tmp/matrix_{how}.dp{filter}.tsv"),
        png = "results_to_upload/{how}/titv_dp{filter}.png",
        rds = temp("tmp/matrix_{how}.dp{filter}.RDS")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 120,
        time_min=lambda wildcards, attempt: attempt * 30,
        tmpdir="tmp"
    params:
        estimate_extra = "nTry=3",
        h=2048,
        w=1536
    wrapper:
        "bio/maftools/trinucleotidematrix_hg38"


rule maftools_oncoplot:
    input:
        rds = "tmp/{how}.dp{filter}.RDS"
    output:
        png = "results_to_upload/{how}/oncoplot_dp{filter}.png"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 5,
        time_min=lambda wildcards, attempt: attempt * 15,
        tmpdir="tmp"
    params:
        non_synonymous = non_synonymous,
        color_named_vector = color_named_vector,
        png_extra="height=1536, width=1536, units='px'",
        oncoplot_extra="top=40"
    wrapper:
        "bio/maftools/oncoplot"


rule maftools_plotmafsummary:
    input:
        rds = "tmp/{how}.dp{filter}.RDS"
    output:
        png = "results_to_upload/{how}/summary_dp{filter}.png"
    params:
        non_synonymous = non_synonymous,
        color_named_vector = color_named_vector,
        png_extra="height=1536, width=1536, units='px'",
        maftools_extra="top=20"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 5,
        time_min=lambda wildcards, attempt: attempt * 15,
        tmpdir="tmp"
    wrapper:
        "bio/maftools/plotmafsummary"


rule load_maf_to_R:
    input:
        maf="tmp/{how}_{filter}.tsv"
    output:
        rds=temp("tmp/{how}.dp{filter}.RDS"),
        summary = temp("tmp/{how}.dp{filter}_summary.txt")
    params:
        summary_prefix="tmp/{how}.dp{filter}"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 5,
        time_min=lambda wildcards, attempt: attempt * 30,
        tmpdir="tmp"
    log:
        "logs/maftools/readmaf/{how}.{filter}.log"
    wrapper:
        "bio/maftools/readmaf"


rule cat_mafs:
    input:
        expand(
            "results_to_upload/{how}/dp{filter}/{sample}_dp{filter}.tsv",
            sample=sample_names,
            allow_missing=True
        )
    output:
        temp("tmp/{how}_{filter}.tsv")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
        time_min=lambda wildcards, attempt: attempt * 15,
        tmpdir="tmp"
    log:
        "logs/maftools/cat/{how}.{filter}.log"
    params:
        sed="'1d'",
        head="-n1"
    shell:
        """head {params.head} {input[0]} > {output} 2> {log} && """
        """for FILEPATH in {input}; do sed {params.sed} "${{FILEPATH}}"; done """
        """>> {output} 2>> {log}"""

##################
### Occurences ###
##################

"""
Count variant occurence
"""
rule variant_occurence_annotate:
    input:
        calls = ["maf/canonical/{sample}.vcf"],
        occurence = "maf/occurences.txt"
    output:
        calls = ["maf/occurence_annotated/{sample}.vcf"]
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
        time_min=lambda wildcards, attempt: attempt * 15,
        tmpdir="tmp"
    log:
        "logs/variant_occurence/uncompress/{sample}.log"
    wrapper:
        "bio/variantoccurence/annotate"


rule concatenate_per_chr_information:
    input:
        expand("maf/{chr}/occurence.txt", chr=[*map(str, range(23)), *range(23), "X", "Y", "MT"])
    output:
        "maf/occurences.txt"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
        time_min=lambda wildcards, attempt: attempt * 15,
        tmpdir="tmp"
    log:
        "logs/variant_occurence/all.log"
    shell:
        "for i in {input}; do sed '1d' ${{i}}; done > {output} 2> {log}"


rule variant_occurence_per_chr:
    input:
        calls=expand(
            "maf/canonical/{sample}.vcf.gz",
            sample=sample_names
        )
    output:
        txt="maf/{chr}/occurence.txt"
    threads: 7
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
        time_min=lambda wildcards, attempt: attempt * 45,
        tmpdir="tmp"
    log:
        "logs/variant_occurence/{chr}.log"
    wrapper:
        "bio/variantoccurence/chromosomes"

##################
### VCF to maf ###
##################

rule gatk_variant_filtration:
    input:
        vcf="maf/occurence_annotated/{sample}.vcf.gz",
        vcf_tbi="maf/occurence_annotated/{sample}.vcf.gz.tbi",
        ref=config["ref"]["fasta"]
    output:
        vcf="results_to_upload/filtered_vcf/{sample}.vcf.gz",
        vcf_idx="results_to_upload/filtered_vcf/{sample}.vcf.gz.tbi"
    message:
        "Filtering VCF for {wildcards.sample}"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 10240,
        time_min=lambda wildcards, attempt: attempt * 25,
        tmpdir="tmp"
    log:
        "logs/gatk/variant_filtration/{sample}.log"
    params:
        filters=lambda wildcards: {
            "AFBelow5pct": "AF < 0.05",
            "DepthBelow5X": "DP < 5",
            "DepthBelow10X": "DP < 10",
            "DepthBelow40X": "DP < 40",
            "AFBelow10pct": "AF < 0.1",
            "KaviarAFBelow5pct": "Kaviar_AF < 0.05",
            "GnomadAFBelow5pct": "dbNSFP_gnomAD_exomes_AF < 0.05",
            "TumorAFBelow5pct": f"FORMAT_{wildcards.sample}_tumor_AF < 0.05",
            #"VarOccAbove105": "VarOcc > 105"
            #"VarOccBelow5": "VarOcc < 5",
        },
        extra="--create-output-variant-md5 --create-output-variant-index "
    wrapper:
        "bio/gatk/variantfiltration"


rule gunzip:
    input:
        "results_to_upload/filtered_vcf/{sample}.vcf.gz"
    output:
        temp("tmp/filtered_vcf/{sample}.vcf")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
        time_min=lambda wildcards, attempt: attempt * 15,
        tmpdir="tmp"
    log:
        "logs/gunzip/{sample}.log"
    params:
        "-c"
    shell:
        "gunzip {params} {input} > {output} 2> {log}"


rule cgc_annotate:
    input:
        vcf="tmp/filtered_vcf/{sample}.vcf",
        cgc="/mnt/beegfs/userdata/t_dayris/Census_allTue_Aug_31_15_11_39_2021.csv"
    output:
        vcf=temp("tmp/cgc/{sample}.vcf")
    threads: 2
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 10240,
        time_min=lambda wildcards, attempt: attempt * 25,
        tmpdir="tmp"
    log:
        "logs/cgc/{sample}.log"
    wrapper:
        "bio/BiGR/cancer_gene_census_annotate"


rule format_to_info:
    input:
        call = "tmp/cgc/{sample}.vcf"
    output:
        call = temp("tmp/f2i/{sample}.vcf")
    message:
        "Moving format fields to info for {wildcards.sample}"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 2048,
        time_min=lambda wildcards, attempt: attempt * 45,
        tmpdir="tmp"
    params:
        normal_sample=lambda wildcards: f"{wildcards.sample}_normal",
        tumor_sample=lambda wildcards: f"{wildcards.sample}_tumor"
    log:
        "logs/vcf_format_to_info/{sample}.log"
    wrapper:
        "bio/BiGR/vcf_format_to_info"


rule oncokb_annotate:
    input:
        vcf="tmp/f2i/{sample}.vcf",
        oncokb="/mnt/beegfs/userdata/t_dayris/OncoKB.csv"
    output:
        vcf="results_to_upload/calling/{sample}.vcf"
    threads: 2
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 10240,
        time_min=lambda wildcards, attempt: attempt * 25,
        tmpdir="tmp"
    log:
        "logs/oncokb/{sample}.log"
    wrapper:
        "bio/BiGR/oncokb_annotate"


rule vcf_to_tsv:
    input:
        call="results_to_upload/calling/{sample}.vcf",
    output:
        tsv=temp("tmp/{sample}.tsv")
    threads: 2
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 10240,
        time_min=lambda wildcards, attempt: attempt * 25,
        tmpdir="tmp"
    log:
        "logs/snpsift/extract_all_fields/{sample}.log"
    params:
        extra="-e '.' -s $'\\t'"
    wrapper:
        "bio/snpsift/extractAllFields"


rule rename_cols:
    input:
        tsv="tmp/{sample}.tsv"
    output:
        tsv=temp("tmp/readable_{sample}.tsv")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 10,
        time_min=lambda wildcards, attempt: attempt * 15,
        tmpdir="tmp"
    log:
        "logs/bigr/rename_cols/{sample}.log"
    params:
        Tumor_Sample_Barcode=lambda wildcards: f"{wildcards.sample}_tumor",
        Matched_Norm_Sample_Barcode=lambda wildcards: f"{wildcards.sample}_normal",
        add_cols=True
    wrapper:
        "bio/BiGR/rename_snpsift_maf_cols"


rule filter_pandas:
    input:
        table = "tmp/translated_maftools_{sample}.tsv" # "tmp/readable_{sample}.tsv"
    output:
        table = "results_to_upload/all/dp{filter}/{sample}_dp{filter}.tsv" # temp("tmp/readable_filtered{filter}_{sample}.tsv")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 10,
        time_min=lambda wildcards, attempt: attempt * 15,
        tmpdir="tmp"
    log:
        "logs/maf/cols/{sample}_{filter}.log"
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
            "Mutect2_Read_depth": "float",
            f"{wildcards.sample}_tumor_AF": "float",
            f"{wildcards.sample}_normal_AF": "float",
            f"{wildcards.sample}_normal_AD_allele2": "float",
            f"{wildcards.sample}_tumor_DP": "float",
            f"{wildcards.sample}_tumor_AD_allele2": "float",
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
            [f"{wildcards.sample}_normal_AD_allele2", "<=", float(wildcards.filter)],
            [f"{wildcards.sample}_tumor_AD_allele2", ">=", float(wildcards.filter)],
            [f"{wildcards.sample}_tumor_DP", ">=", float(wildcards.filter)],
            #["dbNSFP_ExAC_AlleleFrequency", "<=", 0.05],
            ["VarOcc", "<=", len(sample_names) - 1]
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
            ["Filter", "VarOccAbove105"],
            #["Filter", f"DepthBelow{wildcards.filter}X"],
            #["Filter", "AFBelow40pct"]
        ],
        drop_duplicated_lines=True
    wrapper:
        "bio/pandas/filter_table"


rule filter_census:
    input:
        table="results_to_upload/all/dp{filter}/{sample}_dp{filter}.tsv",
    output:
        table="results_to_upload/census_only/dp{filter}/{sample}_dp{filter}.tsv"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 10,
        time_min=lambda wildcards, attempt: attempt * 15,
        tmpdir="tmp"
    log:
        "logs/census_only/{sample}.{filter}.log"
    params:
        contains = [["Filter", "ExistsInCancerGeneCensus"]]
    wrapper:
        "bio/pandas/filter_table"


rule filter_oncokb:
    input:
        table="results_to_upload/all/dp{filter}/{sample}_dp{filter}.tsv",
    output:
        table="results_to_upload/oncokb_only/dp{filter}/{sample}_dp{filter}.tsv"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 10,
        time_min=lambda wildcards, attempt: attempt * 15,
        tmpdir="tmp"
    log:
        "logs/census_only/{sample}.{filter}.log"
    params:
        contains = [["Filter", "ExistsInOncoKB"]]
    wrapper:
        "bio/pandas/filter_table"


rule maftools_translate_variant_classification:
    input:
        maf = "tmp/readable_{sample}.tsv" # "tmp/readable_filtered{filter}_{sample}.tsv"
    output:
        maf = temp("tmp/translated_maftools_{sample}.tsv") #"results_to_upload/all/{sample}_dp{filter}.tsv"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 8,
        time_min=lambda wildcards, attempt: attempt * 15,
        tmpdir="tmp"
    log:
        "logs/maftools/translation/{sample}.log"
    wrapper:
        "bio/maftools/translate_snpeff_maf_effect"


####################
### Compress VCF ###
####################

module compress_index_vcf_meta:
    snakefile: "/mnt/beegfs/pipelines/snakemake-wrappers/meta/bio/compress_index_vcf/test/Snakefile"
    config: config

use rule * from compress_index_vcf_meta
