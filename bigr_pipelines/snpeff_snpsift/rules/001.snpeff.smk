# Annotate with snpeff
rule snpeff:
    input:
        calls="data_input/calls/{sample}.vcf.gz",
        calls_index="data_input/calls/{sample}.vcf.gz.tbi",
        db=config["ref"]["snpeff"]
    output:
        calls=temp("snpeff/calls/{sample}.vcf"),
        stats="snpeff/report/{sample}.html",
        csvstats=temp("snpeff/csvstats/{sample}.csv")
    message: "Annotating {wildcards.sample} with SnpEff"
    threads: 3
    resources:
        mem_mb=lambda wildcard, attempt: attempt * 1024 * 8,
        time_min=lambda wildcard, attempt: attempt * 90,
        tmpdir="tmp"
    params:
        extra=config["snpeff_snpsift"].get("snpeff_extra", "-nodownload -noLog")
    log:
        "logs/snpeff/annotate/{sample}.log"
    wrapper:
        "bio/snpeff/annotate"