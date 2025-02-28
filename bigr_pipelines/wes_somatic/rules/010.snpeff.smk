rule snpeff:
    input:
        calls="mutect2/corrected/{sample}.vcf.gz",
        calls_index="mutect2/corrected/{sample}.vcf.gz.tbi",
        db=config["reference"]["snpeff"],
    output:
        calls=temp("snpeff/calls/{sample}.vcf"),
        stats="snpeff/report/{sample}.html",
        csvstats=temp("snpeff/csvstats/{sample}.csv"),
    threads: 3
    resources:
        mem_mb=get_8gb_per_attempt,
        time_min=get_90min_per_attempt,
        tmpdir=tmp,
    retries: 1
    params:
        extra=config["snpeff"].get("extra", "-nodownload -noLog"),
    log:
        "logs/snpeff/annotate/{sample}.log",
    wrapper:
        "bio/snpeff/annotate"


###############################################################
### Adding a filtering step here, to make the rest of the   ###
### annotation process faster. The filters here shall not   ###
### include annotation-related information                  ###
###                                                         ###
### This is done *after* snpeff in order to have gene names ###
### alongside with reason why they are not kept. This is    ###
### for PI negociations                                     ###
###############################################################


rule gatk_hard_filtering:
    input:
        vcf="snpeff/calls/{sample}.vcf.gz",
        vcf_tbi="snpeff/calls/{sample}.vcf.gz.tbi",
        ref=config["reference"]["fasta"],
        ref_index=config["reference"]["fasta_index"],
        ref_dict=config["reference"]["fasta_dict"],
    output:
        vcf=temp("gatk/variantfiltration/{sample}.vcf.gz"),
        vcf_tbi=temp("gatk/variantfiltration/{sample}.vcf.gz.tbi"),
    threads: 1
    resources:
        mem_mb=get_8gb_per_attempt,
        time_min=get_90min_per_attempt,
        tmpdir=tmp,
    params:
        filters=config["gatk"].get(
            "gatk_filters_quality",
            {
                "DepthBelow10X": "DP < 10",
                "BelowQualByDepth": "QD <= 2.0",
                # "BelowBaseQuality": "QUAL < 30.0",
                "AboveFisherStrandBias": "FS > 60.0",
                "AboveStrandOddsRatio": "SOR > 3.0",
                "BelowMappingQuality": "MQ < 35.0",
                "BelowMQRankSum": "MQRankSum < -12.5",
                "BelowReadPosRankSum": "ReadPosRankSum < -8.0",
            },
        ),
        extra="--create-output-variant-index  --seconds-between-progress-updates 30 --missing-values-evaluate-as-failing false",
    log:
        "logs/gatk/variantfiltration/{sample}.log",
    wrapper:
        "bio/gatk/variantfiltration"


rule bcftools_select_variants_preannot:
    input:
        vcf="gatk/variantfiltration/{sample}.vcf.gz",
        vcf_tbi="gatk/variantfiltration/{sample}.vcf.gz.tbi",
        ref=config["reference"]["fasta"],
        ref_index=config["reference"]["fasta_index"],
        ref_dict=config["reference"]["fasta_dict"],
        regions=config["reference"]["capture_kit_bed"],
    output:
        vcf=temp("bcftools/filter/{sample}.preannot.vcf"),
    threads: 1
    resources:
        mem_mb=get_8gb_per_attempt,
        time_min=get_90min_per_attempt,
        tmpdir=tmp,
    params:
        extra=""# ' --include \'FILTER=="PASS" || FILTER=="."\'',
    log:
        "logs/bcftools/filter/{sample}.pre.annotation.log",
    wrapper:
        "bio/bcftools/filter"
