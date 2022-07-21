gatk_mutect2_somatic_config = {
    "genome": config["ref"]["fasta"],
    "known": config["ref"]["af_only"],
    "bed": config["ref"]["capture_kit_bed"],
    "dbsnp": config["ref"]["dbsnp"],
    "sample_list": design["Sample_id"].to_list(),
    "chrom": config["params"]["chr"]
}


module gatk_mutect2_somatic_meta:
    snakefile: str(worflow_source_dir / ".." / ".." / "meta" / "bio" / "mutect2_somatic" / "test" / "Snakefile")
    config: gatk_mutect2_somatic_config


use rule * from gatk_mutect2_somatic_meta


rule correct_mutect2_vcf:
    input:
        "bcftools/mutect2/{sample}.vcf.gz"
    output:
        temp("mutect2/corrected/{sample}.vcf")
    message:
        "Fixing AS_FilterStrand format error"
        " on {wildcards.sample}"
    threads: 2
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 256,
        time_min=lambda wildcards, attempt: attempt * 20,
        tmpdir="tmp"
    log:
        "logs/mutect2/correct_fields/{sample}.log"
    params:
        fix_as_filterstatus="'s/ID=AS_FilterStatus,Number=A/ID=AS_FilterStatus,Number=1/g'"
    shell:
        "(gunzip -c {input} | "
        "sed {params.fix_as_filterstatus}) "
        "> {output} 2> {log}"