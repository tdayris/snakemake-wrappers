rule extractfields:
    input:
        call="snpsift/fixed/{sample}.vcf.gz",
        call_index=get_tbi("snpsift/fixed/{sample}.vcf.gz"),
    output:
        tsv="snpsift/extractFields/{sample}.tsv",
    message:
        "Making {wildcards.sample} annotated VCF readable"
    threads: 2
    resources:
        mem_mb=lambda wildcards, attempt: min(attempt * 4096, 15360),
        time_min=lambda wildcards, attempt: attempt * 20,
    log:
        "logs/snpsift/extractAllFields/{sample}.log",
    params:
        extra="-s ';' -e '.'",
    wrapper:
        "bio/snpsift/extractAllFields"


rule fix_vcf:
    input:
        vcf=last_vcf,
    output:
        vcf=temp("snpsift/fixed/{sample}.vcf"),
    message:
        "Removing empty fields, trailing ';' and non-canonical chromosomes "
        "for {wildcards.sample}"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
        time_min=lambda wildcards, attempt: attempt * 15,
        tmpdir="tmp",
    log:
        "logs/bigr_scripts/fix_vcf/{sample}.log",
    params:
        default_chr=config["params"]["chr"],
        remove_non_conventional_chromosomes=True,
    wrapper:
        "bio/BiGR/fix_vcf"
