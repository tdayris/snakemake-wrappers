rule snpeff:
    input:
        calls="bcftools/mutect2/{sample}.vcf.gz",
        calls_index="bcftools/mutect2/{sample}.vcf.gz.tbi",
        db=config["reference"]["snpeff"],
    output:
        calls=temp("snpeff/calls/{sample}.vcf"),
        stats="snpeff/report/{sample}.html",
        csvstats=temp("snpeff/csvstats/{sample}.csv"),
    threads: 3
    resources:
        mem_mb=get_8gb_per_attempt,
        time_min=get_90min_per_attempt,
        tmpdir="tmp",
    params:
        extra=config["snpeff"].get("extra", "-nodownload -noLog"),
    log:
        "logs/snpeff/annotate/{sample}.log",
    wrapper:
        "bio/snpeff/annotate"
