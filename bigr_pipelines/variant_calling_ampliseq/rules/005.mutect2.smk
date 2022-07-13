# Fix mutect2 odd Format names
# AS_FilterStatus: Number=1 and not Number=A which violates VCF format
# AD becomes ADM: AD is reserved for Allele Depth, Mutect2 stores
#                 multiple information under "AD" field. 

rule correct_mutect2_vcf:
    input:
        "mutect2/filter_reheaded/{sample}.vcf.gz"
    output:
        temp("mutect2/corrected/{sample}.vcf")
    message:
        "Renaming reserved AD field and fixing AS_FilterStrand format error"
        " on {wildcards.sample}"
    threads: 3
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 256,
        time_min=lambda wildcards, attempt: attempt * 20
    log:
        "logs/mutect2/correct_fields/{sample}.log"
    params:
        rename_ad="'s/=AD;/=ADM;/g'",
        rename_ad_format="'s/:AD:/:ADM:/g'",
        fix_as_filterstatus="'s/ID=AS_FilterStatus,Number=A/ID=AS_FilterStatus,Number=1/g'"
    shell:
        "(gunzip -c {input} | "
        "sed {params.rename_ad} | "
        "sed {params.rename_ad_format} | "
        "sed {params.fix_as_filterstatus}) "
        "> {output} 2> {log}"



# Call Mutect2 meta-wrapper
gatk_mutect2_germline_meta_config = {
    "genome": config["ref"]["fasta"], 
    "known": config["ref"]["af_only"], 
    "bed": config["ref"]["capture_kit_bed"], 
    "dbsnp": config["ref"]["dbsnp"]
}

module gatk_mutect2_germline_meta:
    snakefile: "../../../meta/bio/mutect2_germline/test/Snakefile"
    config: gatk_mutect2_germline_meta_config


use rule * from gatk_mutect2_germline_meta as gatk_mutect2_germline_*
