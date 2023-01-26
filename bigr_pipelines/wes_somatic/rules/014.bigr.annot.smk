#############################
### GLeaves compatibility ###
#############################


###################################
### OncoKB and CancerGeneCensus ###
### Custom annotations          ###
###################################


rule cancer_gene_census_annotate:
    input:
        vcf="bigr/oncokb/{sample}.vcf",
        cgc=config["reference"]["cancer_census"],
    output:
        vcf=temp("bigr/cancer_gene_census/{sample}.vcf"),
    threads: 1
    resources:
        mem_mb=get_6gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir=tmp,
    retries: 1
    log:
        "logs/bigr/cancer_gene_census_annotate/{sample}.log",
    wrapper:
        "bio/BiGR/cancer_gene_census_annotate"


rule oncokb_annotate:
    input:
        vcf="bigr/format_to_info/{sample}.vcf",
        oncokb=config["reference"]["oncokb"],
    output:
        vcf=temp("bigr/oncokb/{sample}.vcf"),
    threads: 1
    resources:
        mem_mb=get_6gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir=tmp,
    retries: 1
    log:
        "logs/bigr/oncokb/{sample}.log",
    wrapper:
        "bio/BiGR/oncokb_annotate"


####################
## Format to info ##
####################


rule format_to_info:
    input:
        call="vcftools/mane/{sample}.vcf.gz",
    output:
        call=temp("bigr/format_to_info/{sample}.vcf"),
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_15min_per_attempt,
        tmpdir=tmp,
    retries: 1
    log:
        "logs/bigr/format_to_info/{sample}.log",
    params:
        extra=config["bigr_additionals"].get("format_to_info", ""),
    wrapper:
        "bio/BiGR/vcf_format_to_info"
