rule filter_tsv:
    input:
        "snpsift/extractFields/{sample}.tsv",
    output:
        table=protected("data_output/TSV/{sample}.tsv"),
        xlsx=protected("data_output/XLSX/{sample}.xlsx"),
    threads: 1
    resources:
        mem_mb=get_8gb_per_attempt,
        time_min=get_20min_per_attempt,
        tmpdir="tmp",
    retries: 1
    params:
        drop_duplicated_lines=True,
        contains=[["Filter", "PASS"]],
        prefixes=[["Chromosome", "chr"]],
        new_cols=lambda wildcards: [
            ["Mutect2_Allele_Frequency", "=", f"{wildcards.sample}_tumor_AF"]
        ],
        keep_column=lambda wildcards: config["table_cols"]
        + [
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
    wrapper:
        "bio/pandas/filter_table"


rule extractfields:
    input:
        call="data_output/VCF/{sample}.vcf.gz",
        call_index="data_output/VCF/{sample}.vcf.gz.tbi",
    output:
        tsv="snpsift/extractFields/{sample}.tsv",
    threads: 2
    resources:
        mem_mb=get_6gb_per_attempt,
        time_min=get_20min_per_attempt,
        tmpdir="tmp",
    retries: 1
    log:
        "logs/snpsift/extractAllFields/{sample}.log",
    params:
        extra=config["snpsift"].get("extract_all_fields", "-s ';' -e '.'"),
    wrapper:
        "bio/snpsift/extractAllFields"


rule fix_vcf:
    input:
        vcf=last_vcf,
    output:
        vcf=temp("snpsift/fixed/{sample}.vcf"),
    threads: 1
    resources:
        mem_mb=get_10gb_per_attempt,
        time_min=get_15min_per_attempt,
        tmpdir="tmp",
    retries: 1
    log:
        "logs/bigr_scripts/fix_vcf/{sample}.log",
    params:
        default_chr=config["reference"]["chr"],
        remove_non_conventional_chromosomes=config["bigr_additionals"].get(
            "remove_non_conventional_chromosomes",
            True,
        ),
    wrapper:
        "bio/BiGR/fix_vcf"


rule gleaves_compatibility:
    input:
        "snpsift/fixed/{sample}.vcf",
    output:
        vcf=temp("gleaves/corrected/{sample}.vcf.gz"),
        vcf_tbi=temp("gleaves/corrected/{sample}.vcf.gz.tbi"),
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir="tmp",
    retries: 1
    log:
        "logs/BiGR/gleaves_compatibility/{sample}.log",
    wrapper:
        "bio/BiGR/gleaves_compatibility"


rule gatk_variant_filtration:
    input:
        vcf="gleaves/corrected/{sample}.vcf.gz",
        vcf_index="gleaves/corrected/{sample}.vcf.gz.tbi",
        ref=config["reference"]["fasta"],
        ref_idx=config["reference"]["fasta_index"],
        ref_dict=config["reference"]["fasta_dict"],
    output:
        vcf=protected("data_output/VCF/{sample}.vcf.gz"),
        vcf_tbi=protected("data_output/VCF/{sample}.vcf.gz.tbi"),
    threads: 1
    resources:
        mem_mb=get_8gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir="tmp",
    retries: 1
    log:
        "logs/gatk/variant_filtration/{sample}.log",
    params:
        filters=config["gatk"].get("filters", {"DepthBelow10X": "DP < 10"}),
        extra="--create-output-variant-index --create-output-variant-md5",
    wrapper:
        "bio/gatk/variantfiltration"