rule cnv_facets:
    input:
        tumor_bam="sambamba/markdup/{sample}_tumor.bam",
        normal_bam="sambamba/markdup/{sample}_normal.bam",
        vcf=config["reference"]["dbsnp"],
        vcf_index=config["reference"]["dbsnp_tbi"],
        bed=config["reference"]["capture_kit_bed"],
    output:
        vcf=temp("facets/{sample}/{sample}.vcf.gz"),
        profile=temp("facets/{sample}/{sample}.cnv.png"),
        coverage=temp("facets/{sample}/{sample}.cov.pdf"),
        spider=temp("facets/{sample}/{sample}.spider.pdf"),
        pileup=temp("facets/{sample}/{sample}.csv.gz"),
    threads: min(config.get("max_threads", 20), 10)
    resources:
        mem_mb=get_30gb_and_10gb_per_attempt,
        time_min=get_2h_per_attempt,
        tmpdir=tmp,
    retries: 1
    params:
        extra=config["cnv_facets"].get(
            "extra", "--snp-count-orphans --gbuild hg38 --nbhd-snp 250"
        ),
        prefix="facets/{sample}/{sample}",
    log:
        "logs/facets/cnv/{sample}.log",
    wrapper:
        str(wrapper_prefix / "bio" / "facets" / "cnv")
