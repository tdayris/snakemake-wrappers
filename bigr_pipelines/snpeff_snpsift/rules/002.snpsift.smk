
rule snpsift_exac:
    input:
        call = "snpsift/gnomad/{sample}.vcf",
        database = config["ref"]["exac"],
        database_idx = config["ref"]["exac"] + ".tbi"
    output:
        call = temp("snpsift/exac/{sample}.vcf")
    message:
        "Annotating {wildcards.sample} with ExAC"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1020 + 4096,
        time_min=lambda wildcards, attempt: attempt * 75,
        tmpdir="tmp"
    params:
        extra=config["snpeff_snpsift"].get(
            "snpsift_exac",
            "-name 'gnomAD_exomes_' -exists 'ExistsInExac' -tabix -noDownload -noLog"
        )
    log:
        "logs/snpsift/exac/{sample}.log"
    wrapper:
        "bio/snpsift/annotate"



rule snpsift_gnomad:
    input:
        call = "snpsift/clinvar/{sample}.vcf",
        database = config["ref"]["gnomad"],
        database_idx = config["ref"]["gnomad"] + ".tbi"
    output:
        call = temp("snpsift/gnomad/{sample}.vcf")
    message:
        "Annotating {wildcards.sample} with gnomAD"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1020 + 4096,
        time_min=lambda wildcards, attempt: attempt * 45,
        tmpdir="tmp"
    params:
        extra=config["snpeff_snpsift"].get(
            "snpsift_gnomad",
            "-name 'gnomad_' -exists 'ExistsInGnomAD' -tabix -noDownload -noLog"
        )
    log:
        "logs/snpsift/gnomad/{sample}.log"
    wrapper:
        "bio/snpsift/annotate"


rule snpsift_clinvar:
    input:
        call = "snpsift/dbnsfp/{sample}.vcf",
        database = config["ref"]["clinvar"],
        database_idx = config["ref"]["clinvar"] + ".tbi"
    output:
        call = temp("snpsift/clinvar/{sample}.vcf")
    message:
        "Annotating {wildcards.sample} with ClinVar"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1020 + 4096,
        time_min=lambda wildcards, attempt: attempt * 45,
        tmpdir="tmp"
    params:
        extra=config["snpeff_snpsift"].get(
            "snpsift_clinvar",
            "-name 'clinvar_' -exists 'ExistsInClinVar' -tabix -noDownload -noLog"
        )
    log:
        "logs/snpsift/clinvar/{sample}.log"
    wrapper:
        "bio/snpsift/annotate"


rule snpsift_gwascat:
    input:
        call = "snpsift/dbnsfp/{sample}.vcf",
        gwascat = config["ref"]["gwascat"]
    output:
        call = temp("snpsift/gwascat/{sample}.vcf")
    message:
        "Annotating {wildcards.sample} with GWAS Catalog"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1020 + 4096,
        time_min=lambda wildcards, attempt: attempt * 45,
        tmpdir="tmp"
    params:
        extra = config["snpeff_snpsift"].get(
            "snpsift_gwascat", "-noDownload -noLog"
        )
    log:
        "logs/snpsift/gwascat/{sample}.log"
    wrapper:
        "bio/snpsift/gwascat"


rule snpsift_dbnsfp:
    input:
        call = "snpsift/cosmic/{sample}.vcf",
        dbNSFP = config["ref"]["dbnsfp"],
        dbNSFP_tbi = config["ref"]["dbnsfp"] + ".tbi"
    output:
        call = temp("snpsift/dbnsfp/{sample}.vcf")
    message:
        "Annotating {wildcards.sample} with dbNSFP"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 20480 + 10240,
        time_min=lambda wildcards, attempt: attempt * 45,
        tmpdir="tmp"
    params:
        extra=config["snpeff_snpsift"].get(
            "snpsift_dbnsfp",
            "-noDownload -noLog -n -f 'hg18_chr,hg18_pos(1-based)'"
        )
    log:
        "logs/snpsift/dbnsfp/{sample}.log"
    wrapper:
        "bio/snpsift/dbnsfp"


rule snpsift_cosmic:
    input:
        call="snpsift/kaviar/{sample}.vcf",
        database=config["ref"]["cosmic"],
        database_idx=config["ref"]["cosmic"] + ".tbi"
    output:
        call=temp("snpsift/cosmic/{sample}.vcf")
    message:
        "Annotating {wildcards.sample} with COSMIC"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1020 + 4096,
        time_min=lambda wildcards, attempt: attempt * 120,
        tmpdir="tmp"
    params:
        extra=config["snpeff_snpsift"].get(
            "snpsift_cosmic",
            "-name 'cosmic_' -exists 'ExistsInCosmic' -tabix -noDownload -noLog"
        )
    log:
        "logs/snpsift/cosmic/{sample}.log"
    wrapper:
        "bio/snpsift/annotate"


rule snpsift_kaviar:
    input:
        call="snpsift/gmt/{sample}.vcf",
        database=config["ref"]["kaviar"],
        database_idx=config["ref"]["kaviar"] + ".tbi",
    output:
        call=temp("snpsift/kaviar/{sample}.vcf")
    message:
        "Annotating {wildcards.sample} with Kaviar"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1020 + 4096,
        time_min=lambda wildcards, attempt: attempt * 45,
        tmpdir="tmp"
    log:
        "logs/snpsift/kaviar/{sample}.log"
    params:
        extra=config["snpeff_snpsift"].get(
            "snpsift_kaviar",
            "-name 'Kaviar_' -exists 'ExistsInKaviar' -tabix -noDownload -noLog"
        )
    wrapper:
        "bio/snpsift/annotate"


rule snpsift_gmt:
    input:
        call = "snpsift/dbsnp/{sample}.vcf",
        gmt = config["ref"]["gmt"]
    output:
        call = temp("snpsift/gmt/{sample}.vcf")
    message:
        "Annotating {wildcards.sample} with MSigDB"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1020 + 4096,
        time_min=lambda wildcards, attempt: attempt * 45,
        tmpdir="tmp"
    log:
        "logs/snpsift/gmt/{sample}.log"
    wrapper:
        "bio/snpsift/genesets"


rule snpsift_dbsnp:
    input:
        call="snpsift/vartype/{sample}.vcf",
        database=config["ref"]["dbsnp"]
    output:
        call=temp("snpsift/dbsnp/{sample}.vcf")

    message:
        "Annotating {wildcards.sample} with dbSNP"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1020 + 4096,
        time_min=lambda wildcards, attempt: attempt * 45,
        tmpdir="tmp"
    log:
        "logs/snpsift/dbsnp/{sample}.log"
    params:
        extra=config["snpeff_snpsift"].get(
            "snpsift_dbsnp",
            "-name 'dbSNP_' -exists 'ExistsInDBsnp' -tabix -noDownload -noLog"
        )
    wrapper:
        "bio/snpsift/annotate"


rule snpsift_vartype:
    input:
        call="snpeff/calls/{sample}.vcf",
        #vcf_tbi="snpeff/calls/{sample}.vcf.gz.tbi"
    output:
        call=temp("snpsift/vartype/{sample}.vcf")
    message:
        "Annotating variant types in {wildcards.sample}"
    threads: 2
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1020 + 4096,
        time_min=lambda wildcards, attempt: attempt * 45,
        tmpdir="tmp"
    log:
        "logs/snpsift/varType/{sample}.log"
    wrapper:
        "bio/snpsift/varType"