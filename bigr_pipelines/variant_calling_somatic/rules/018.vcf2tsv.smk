rule extractfields:
    input:
        call="snpsift/fixed/{sample}.vcf.gz",
        call_index="snpsift/fixed/{sample}.vcf.gz.tbi",
    output:
        tsv="snpsift/extractFields/{sample}.tsv",
    threads: 2
    resources:
        mem_mb=get_6gb_per_attempt,
        time_min=get_20min_per_attempt,
    log:
        "logs/snpsift/extractAllFields/{sample}.log",
    params:
        extra=config["snpsift"].get("extract_all_fields", "-s ';' -e '.'"),
    wrapper:
        "bio/snpsift/extractAllFields"


rule fix_vcf:
    input:
        vcf=last_vcf,
    output:
        vcf=temp("snpsift/fixed/{sample}.vcf"),
    threads: 1
    resources:
        mem_mb=get_10gb_per_attempt,
        time_min=get_15min_per_attempt,
        tmpdir="tmp",
    log:
        "logs/bigr_scripts/fix_vcf/{sample}.log",
    params:
        default_chr=config["reference"]["chr"],
        remove_non_conventional_chromosomes=config["bigr_additionals"].get(
            "remove_non_conventional_chromosomes"
            True,
        )
    wrapper:
        "bio/BiGR/fix_vcf"
