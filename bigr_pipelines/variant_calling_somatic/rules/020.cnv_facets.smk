rule cnv_facets:
    input:
        #pileup="samtools/mpileup/{sample}.chr.mpileup.gz",
        tumor_bam="sambamba/markdup/{sample}_tumor.bam",
        normal_bam="sambamba/markdup/{sample}_normal.bam",
        vcf=config["ref"]["dbsnp"],
        vcf_index=get_tbi(config["ref"]["dbsnp"]),
        bed=config["ref"]["capture_kit_bed"],
    output:
        vcf="facets/{sample}/{sample}.vcf.gz",
        profile="facets/{sample}/{sample}.cnv.png",
        coverage="facets/{sample}/{sample}.cov.pdf",
        spider="facets/{sample}/{sample}.spider.pdf",
        pileup="facets/{sample}/{sample}.csv.gz",
    message:
        "Searching for CNV in {wildcards.sample} with Facets"
    threads: 10
    resources:
        mem_mb=lambda wildcards, attempt: (attempt * 1024 * 20) + 102400,
        time_min=lambda wildcards, attempt: attempt * 60 * 2,
        tmpdir="tmp",
    params:
        extra=config.get(
            "facets_extra", "--snp-count-orphans --gbuild hg38 --nbhd-snp 250"
        ),
        prefix="facets/{sample}/{sample}",
    log:
        "logs/facets/cnv/{sample}.log",
    wrapper:
        "bio/facets/cnv"
