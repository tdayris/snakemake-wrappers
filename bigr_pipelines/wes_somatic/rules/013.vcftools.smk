#####################
## MANE annotation ##
#####################


rule additional_headers_mane:
    output:
        temp("mane/description.txt"),
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_15min_per_attempt,
        tmpdir=tmp,
    conda:
        str(workflow_source_dir / "envs" / "bash.yaml")
    group:
        "additional_headers"
    log:
        "logs/mane/description.log",
    params:
        'key=INFO,ID=MANE_NCBI_GeneID,Number=1,Type=String,Description="NCBI_GeneID field from MANE"',
        'key=INFO,ID=MANE_Ensembl_Gene,Number=1,Type=String,Description="Ensembl_Gene field from MANE"',
        'key=INFO,ID=MANE_HGNC_ID,Number=1,Type=String,Description="HGNC_ID field from MANE"',
        'key=INFO,ID=MANE_symbol,Number=1,Type=String,Description="symbol field from MANE"',
        'key=INFO,ID=MANE_name,Number=1,Type=String,Description="name field from MANE"',
        'key=INFO,ID=MANE_RefSeq_nuc,Number=1,Type=String,Description="RefSeq_nuc field from MANE"',
        'key=INFO,ID=MANE_RefSeq_prot,Number=1,Type=String,Description="RefSeq_prot field from MANE"',
        'key=INFO,ID=MANE_Ensembl_nuc,Number=1,Type=String,Description="Ensembl_nuc field from MANE"',
        'key=INFO,ID=MANE_Ensembl_prot,Number=1,Type=String,Description="Ensembl_prot field from MANE"',
        'key=INFO,ID=MANE_MANE_status,Number=1,Type=String,Description="MANE_status field from MANE"',
        'key=INFO,ID=MANE_GRCh38_chr,Number=1,Type=String,Description="GRCh38_chr field from MANE"',
        'key=INFO,ID=MANE_chr_start,Number=1,Type=String,Description="chr_start field from MANE"',
        'key=INFO,ID=MANE_chr_end,Number=1,Type=String,Description="chr_end field from MANE"',
        'key=INFO,ID=MANE_chr_strand,Number=1,Type=String,Description="chr_strand field from MANE"',
    shell:
        'for PARAM in {params}; do echo "${{PARAM}}"; done > {output} 2> {log}'


rule vcftools_annotate_mane:
    input:
        vcf="vcftools/revel/{sample}.vcf.gz",
        vcf_tbi="vcftools/revel/{sample}.vcf.gz.tbi",
        annotation=config["reference"]["mane"],
        description="mane/description.txt",
    output:
        vcf=temp("vcftools/mane/{sample}.vcf.gz"),
    threads: 4
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_15min_per_attempt,
        tmpdir=tmp,
    retries: 1
    log:
        "logs/vcftools/annotate/{sample}.mane.log",
    params:
        extra=config["vcftools"].get(
            "mane",
        "--columns INFO/MANE_NCBI_GeneID,INFO/MANE_Ensembl_Gene,"
            "INFO/MANE_HGNC_ID,INFO/MANE_symbol,INFO/MANE_name,"
            "INFO/MANE_RefSeq_nuc,INFO/MANE_RefSeq_prot,INFO/MANE_Ensembl_nuc,"
            "INFO/MANE_Ensembl_prot,INFO/MANE_MANE_status,CHROM,FROM,TO,"
            "INFO/MANE_chr_strand",
        ),
    wrapper:
        str(wrapper_prefix / "bio" / "vcftools" / "annotate")


######################
## Revel annotation ##
######################


rule additional_headers_revel:
    output:
        temp("revel/description.txt"),
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_15min_per_attempt,
        tmpdir=tmp,
    conda:
        str(workflow_source_dir / "envs" / "bash.yaml")
    group:
        "additional_headers"
    log:
        "logs/mane/description.log",
    params:
        'key=INFO,ID=REVEL_aaref,Number=1,Type=String,Description="Reference Amino Acid from REVEL"',
        'key=INFO,ID=REVEL_aaalt,Number=1,Type=String,Description="Alternative Amino Acid from REVEL"',
        'key=INFO,ID=REVEL,Number=1,Type=String,Description="REVEL score"',
        'key=INFO,ID=REVEL_Ensembl_transcriptid,Number=1,Type=String,Description="Ensemble transcript id from REVEL"',
    shell:
        'for PARAM in {params}; do echo "${{PARAM}}"; done > {output} 2> {log}'


rule vcftools_annotate_revel:
    input:
        vcf="vcftools/mistic/{sample}.vcf.gz",
        vcf_tbi="vcftools/mistic/{sample}.vcf.gz.tbi",
        annotation=config["reference"]["revel"],
        description="revel/description.txt",
    output:
        vcf=temp("vcftools/revel/{sample}.vcf.gz"),
    threads: 4
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_15min_per_attempt,
        tmpdir=tmp,
    retries: 1
    log:
        "logs/vcftools/annotate/{sample}.revel.log",
    params:
        extra=config["vcftools"].get(
            "revel",
        "--columns CHROM,POS,-,REF,ALT,INFO/REVEL_aaref,INFO/REVEL_aaalt,INFO/REVEL,INFO/REVEL_Ensembl_transcriptid"
            if config["reference"]["ncbi_build"].lower() == "grch38"
            else "--columns CHROM,-,POS,REF,ALT,INFO/AAREF,INFO/AAALT,INFO/REVEL,INFO/REVEL_Ensembl_transcriptid",
        ),
    wrapper:
        str(wrapper_prefix / "bio" / "vcftools" / "annotate")


#######################
## Mistic annotation ##
#######################


rule additional_headers_mistic:
    output:
        temp("mistic/description.txt"),
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_15min_per_attempt,
        tmpdir=tmp,
    conda:
        str(workflow_source_dir / "envs" / "bash.yaml")
    group:
        "additional_headers"
    log:
        "logs/mane/description.log",
    params:
        'key=INFO,ID=MISTIC_score,Number=1,Type=String,Description="MISTIC high sensitivity prediction for exome analysis"',
        'key=INFO,ID=MISTIC_pred,Number=1,Type=String,Description="MISTIC prediction for global performance"',
    shell:
        'for PARAM in {params}; do echo "${{PARAM}}"; done > {output} 2> {log}'


rule vcftools_annotate_mistic:
    input:
        vcf="snpsift/exac/{sample}.vcf.gz",
        vcf_tbi="snpsift/exac/{sample}.vcf.gz.tbi",
        annotation=config["reference"]["mistic"],
        description="mistic/description.txt",
    output:
        vcf=temp("vcftools/mistic/{sample}.vcf.gz"),
    threads: 4
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_15min_per_attempt,
        tmpdir=tmp,
    retries: 1
    log:
        "logs/vcftools/annotate/{sample}.mistic.log",
    params:
        extra=config["vcftools"].get(
            "mistic", "--columns CHROM,POS,REF,ALT,INFO/MISTIC_score,INFO/MISTIC_pred"
        ),
    wrapper:
        str(wrapper_prefix / "bio" / "vcftools" / "annotate")
