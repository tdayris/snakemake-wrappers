###################################
### OncoKB and CancerGeneCensus ###
### Custom annotations          ###
###################################


rule cancer_gene_census_annotate:
    input:
        vcf="bigr/oncokb/{sample}.vcf",
        cgc=config["ref"]["cancer_census"]
    output:
        vcf=temp("bigr/cancer_gene_census/{sample}.vcf")
    message:
        "Adding CancerGeneCensus annotation in {wildcards.sample}"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 5,
        time_min=lambda wildcards, attempt: attempt * 25,
        tmpdir="tmp"
    log:
        "logs/bigr/cancer_gene_census_annotate/{sample}.log"
    wrapper:
        "bio/BiGR/cancer_gene_census_annotate"


rule oncokb_annotate:
    input:
        #vcf="snpsift/clinvar/{sample}.vcf",
        vcf="bigr/format_to_info/{sample}.vcf",
        oncokb=config["ref"]["oncokb"]
    output:
        vcf=temp("bigr/oncokb/{sample}.vcf")
    message:
        "Adding OncoKB annotation in {wildcards.sample}"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 6,
        time_min=lambda wildcards, attempt: attempt * 35,
        tmpdir="tmp"
    log:
        "logs/bigr/oncokb/{sample}.log"
    wrapper:
        "bio/BiGR/oncokb_annotate"


####################
## Format to info ##
####################

rule format_to_info:
    input:
        call="vcftools/mane/{sample}.vcf.gz"
    output:
        call="bigr/format_to_info/{sample}.vcf"
    message:
        "Annotating {wildcards.sample} with clear Format descriptions"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
        time_min=lambda wildcards, attempt: attempt * 15,
        tmpdir="tmp"
    log:
        "logs/bigr/format_to_info/{sample}.log"
    params:
        extra = ""
    wrapper:
        "bio/BiGR/vcf_format_to_info"