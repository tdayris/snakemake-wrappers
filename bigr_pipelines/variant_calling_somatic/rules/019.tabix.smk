"""
This rule indexes the pbgzipped vcf file
"""


rule tabix_index:
    input:
        "{tool}/{subcommand}/{sample}.vcf.gz",
    output:
        "{tool}/{subcommand}/{sample}.vcf.gz.tbi",
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir="tmp",
    retries: 1
    params:
        "-p vcf",
    log:
        "logs/{tool}/{subcommand}/tabix/index/{sample}.log",
    wrapper:
        "bio/tabix"


""" 
This rule compress a provided VCF file with pbgzip
"""


rule pbgzip_compress:
    input:
        "{tool}/{subcommand}/{sample}.vcf",
    output:
        "{tool}/{subcommand}/{sample}.vcf.gz",
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_1h_per_attempt,
        tmpdir="tmp",
    retries: 1
    params:
        "",
    log:
        "logs/{tool}/{subcommand}/pbgzip/{sample}.log",
    wrapper:
        "bio/bcftools/view"
