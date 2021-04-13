from pathlib import Path

rule all:
    input:
        expand(
            "snpeff/gwascat/{sample}.vcf",
            sample=[
                str(f.name)[:-len(".vcf.gz")]
                for f in Path("calls").iterdir() 
                if str(f.name).lower().endswith(".vcf.gz")
            ]
        )

rule test_snpsift_vartype:
    input:
        vcf="calls/{sample}.vcf.gz"
    output:
        vcf=temp("snpsift/vartype/{sample}.vcf")
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1020 + 4096,
        time_min=lambda wildcards, attempt: attempt * 45
    log:
        "logs/snpsift/varType/{sample}.log"
    wrapper:
        "/bio/snpsift/varType"


rule snpeff_annotate:
    input:
        calls="snpsift/vartype/{sample}.vcf",
        db="/mnt/beegfs/database/bioinfo/Index_DB/SnpEff/GRCh37.75/"
    output:
        calls=temp("snpeff/{sample}.vcf"),   # annotated calls (vcf, bcf, or vcf.gz)
        stats="snpeff/{sample}.html",        # summary statistics (in HTML), optional
        csvstats="snpeff/{sample}.csv"       # summary statistics in CSV, optional
    log:
        "logs/snpeff/{sample}.log"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1020 + 4096,
        time_min=lambda wildcards, attempt: attempt * 45
    wrapper:
        "/bio/snpeff/annotate"


rule test_snpsift_gmt:
    input:
        call = "snpeff/{sample}.vcf",
        gmt = "/mnt/beegfs/database/bioinfo/Index_DB/MSigDB/GeneSets.gmt"
    output:
        call = temp("snpsift/gmt/{sample}.vcf")
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1020 + 4096,
        time_min=lambda wildcards, attempt: attempt * 45
    wrapper:
        "/bio/snpsift/genesets"


rule snpsift_kaviar:
    input:
        call="snpeff/{sample}.vcf",
        database="/mnt/beegfs/database/bioinfo/Index_DB/Kaviar/HG19/Kaviar-160204-Public/vcfs/Kaviar-160204-Public-hg19-trim.vcf.gz"
    output:
        call=temp("snpeff/kaviar/{sample}.vcf")
    log:
        "logs/snpsift/kaviar/{sample}.log"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1020 + 4096,
        time_min=lambda wildcards, attempt: attempt * 45
    wrapper:
        "/bio/snpsift/annotate"


rule snpsift_dbsnp:
    input:
        call="snpeff/kaviar/{sample}.vcf",
        database="/mnt/beegfs/database/bioinfo/Index_DB/dbSNP/hg19/144_20150605/All_20150605.vcf.gz"
    output:
        call=temp("snpeff/dbsnp/{sample}.vcf")
    log:
        "logs/snpsift/dbsnp/{sample}.log"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1020 + 4096,
        time_min=lambda wildcards, attempt: attempt * 45
    wrapper:
        "/bio/snpsift/annotate"
        
        
rule snpsift_cosmic:
    input:
        call="snpeff/kaviar/{sample}.vcf",
        database="/mnt/beegfs/database/bioinfo/COSMIC/73_20150629/CosmicCodingMuts.vcf"
    output:
        call=temp("snpeff/cosmic/{sample}.vcf")
    log:
        "logs/snpsift/cosmic/{sample}.log"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1020 + 4096,
        time_min=lambda wildcards, attempt: attempt * 45
    wrapper:
        "/bio/snpsift/annotate"


rule snpsift_gwascat:
    input:
        call = "snpeff/cosmic/{sample}.vcf",
        gwascat = "/mnt/beegfs/database/bioinfo/Index_DB/GWASCatalog/gwas_catalog_v1.0.2-studies_r2020-05-03.tsv"
    output:
        call = "snpeff/gwascat/{sample}.vcf"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1020 + 4096,
        time_min=lambda wildcards, attempt: attempt * 45
    wrapper:
        "/bio/snpsift/gwascat"
