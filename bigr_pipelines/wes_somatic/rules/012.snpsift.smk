rule snpsift_gwascat:
    input:
        call="snpsift/dbvar/{sample}.vcf",
        gwascat=config["reference"]["gwascat"],
    output:
        call=temp("snpsift/gwascat/{sample}.vcf"),
    message:
        "Annotating {wildcards.sample} with GWAS Catalog"
    threads: 1
    resources:
        mem_mb=get_6gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir="tmp",
    retries: 1
    params:
        extra=config.get("snpsift_gwascat", "-noDownload -noLog"),
    log:
        "logs/snpsift/gwascat/{sample}.log",
    wrapper:
        "bio/snpsift/gwascat"


rule snpsift_dbvar:
    input:
        call="snpsift/exac/{sample}.vcf",
        database=config["reference"]["dbvar"],
        database_idx=config["reference"]["dbvar_tbi"],
    output:
        call=temp("snpsift/dbvar/{sample}.vcf"),
    threads: 1
    resources:
        mem_mb=get_20gb_and_10gb_per_attempt,
        time_min=get_5h_per_attempt,
        tmpdir="tmp",
    retries: 1
    params:
        extra=config["snpsift"].get("dbvar", "-name 'dbVar_' -tabix -noDownload -noLog"),
    log:
        "logs/snpsift/dbvar/{sample}.log",
    wrapper:
        "bio/snpsift/annotate"


rule snpsift_exac:
    input:
        call="snpsift/gnomad/{sample}.vcf",
        database=config["reference"]["exac"],
        database_idx=config["reference"]["exac_tbi"],
    output:
        call=temp("snpsift/exac/{sample}.vcf"),
    threads: 1
    resources:
        mem_mb=get_6gb_per_attempt,
        time_min=get_75min_per_attempt,
        tmpdir="tmp",
    retries: 1
    params:
        extra=config["snpsift"].get(
            "exac", "-name 'gnomAD_exomes_' -tabix -noDownload -noLog"
        ),
    log:
        "logs/snpsift/exac/{sample}.log",
    wrapper:
        "bio/snpsift/annotate"


rule snpsift_gnomad:
    input:
        call="snpsift/clinvar/{sample}.vcf",
        database=config["reference"]["gnomad"],
        database_idx=config["reference"]["gnomad_tbi"],
    output:
        call=temp("snpsift/gnomad/{sample}.vcf"),
    threads: 1
    resources:
        mem_mb=get_6gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir="tmp",
    retries: 1
    params:
        extra=config["snpsift"].get(
            "gnomad", "-name 'gnomad_' -tabix -noDownload -noLog"
        ),
    log:
        "logs/snpsift/gnomad/{sample}.log",
    wrapper:
        "bio/snpsift/annotate"


rule snpsift_clinvar:
    input:
        call="snpsift/dbnsfp/{sample}.vcf",
        database=config["reference"]["clinvar"],
        database_idx=config["reference"]["clinvar_tbi"],
    output:
        call=temp("snpsift/clinvar/{sample}.vcf"),
    threads: 1
    resources:
        mem_mb=get_6gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir="tmp",
    retries: 1
    params:
        extra=config["snpsift"].get(
            "clinvar", "-name 'clinvar_' -tabix -noDownload -noLog"
        ),
    log:
        "logs/snpsift/clinvar/{sample}.log",
    wrapper:
        "bio/snpsift/annotate"


rule snpsift_dbnsfp:
    input:
        call="snpsift/cosmic/{sample}.vcf",
        dbNSFP=config["reference"]["dbnsfp"],
    output:
        call=temp("snpsift/dbnsfp/{sample}.vcf"),
    threads: 1
    resources:
        mem_mb=get_20gb_and_10gb_per_attempt,
        time_min=lambda wildcards, attempt: attempt * 45,
        tmpdir="tmp",
    retries: 1
    params:
        extra=config["snpsift"].get(
            "dbnsfp", "-noDownload -noLog -n -f 'hg18_chr,hg18_pos(1-based)'"
        ),
    log:
        "logs/snpsift/dbnsfp/{sample}.log",
    wrapper:
        "bio/snpsift/dbnsfp"


rule snpsift_cosmic:
    input:
        call="snpsift/kaviar/{sample}.vcf",
        database=config["reference"]["cosmic"],
        database_idx=config["reference"]["cosmic_tbi"],
    output:
        call=temp("snpsift/cosmic/{sample}.vcf"),
    threads: 1
    resources:
        mem_mb=get_6gb_per_attempt,
        time_min=get_2h_per_attempt,
        tmpdir="tmp",
    retries: 1
    params:
        extra=config["snpsift"].get(
            "cosmic", "-name 'cosmic_' -tabix -noDownload -noLog"
        ),
    log:
        "logs/snpsift/cosmic/{sample}.log",
    wrapper:
        "bio/snpsift/annotate"


rule snpsift_kaviar:
    input:
        call="snpsift/gmt/{sample}.vcf",
        database=config["reference"]["kaviar"],
        database_idx=config["reference"]["kaviar_tbi"],
    output:
        call=temp("snpsift/kaviar/{sample}.vcf"),
    message:
        "Annotating {wildcards.sample} with Kaviar"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1020 + 4096,
        time_min=lambda wildcards, attempt: attempt * 45,
        tmpdir="tmp",
    retries: 1
    log:
        "logs/snpsift/kaviar/{sample}.log",
    params:
        extra=config["snpsift"].get(
            "kaviar", "-name 'Kaviar_' -tabix -noDownload -noLog"
        ),
    wrapper:
        "bio/snpsift/annotate"


rule snpsift_gmt:
    input:
        call="snpsift/dbsnp/{sample}.vcf",
        gmt=config["reference"]["gmt"],
    output:
        call=temp("snpsift/gmt/{sample}.vcf"),
    threads: 1
    resources:
        mem_mb=get_6gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir="tmp",
    retries: 1
    log:
        "logs/snpsift/gmt/{sample}.log",
    wrapper:
        "bio/snpsift/genesets"


rule snpsift_dbsnp:
    input:
        call="snpsift/vartype/{sample}.vcf",
        database=config["reference"]["dbsnp"],
        database_idx=config["reference"]["dbsnp_tbi"],
    output:
        call=temp("snpsift/dbsnp/{sample}.vcf"),
    threads: 1
    resources:
        mem_mb=get_6gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir="tmp",
    retries: 1
    log:
        "logs/snpsift/dbsnp/{sample}.log",
    params:
        extra=config["snpsift"].get("dbsnp", "-name 'dbSNP_' -tabix -noDownload -noLog"),
    wrapper:
        "bio/snpsift/annotate"


rule snpsift_vartype:
    input:
        call="splice_ai/annot/{sample}.vcf.gz",
        calls_tbi="splice_ai/annot/{sample}.vcf.gz.tbi",
    output:
        call=temp("snpsift/vartype/{sample}.vcf"),
    threads: 1
    resources:
        mem_mb=get_6gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir="tmp",
    retries: 1
    log:
        "logs/snpsift/varType/{sample}.log",
    wrapper:
        "bio/snpsift/varType"
