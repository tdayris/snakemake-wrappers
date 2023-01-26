"""
This rule indexes the pbgzipped vcf file
"""


rule tabix_index:
    input:
        "{tool}/{subcommand}/{sample}.vcf.gz",
    output:
        temp("{tool}/{subcommand}/{sample}.vcf.gz.tbi"),
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir=tmp,
    retries: 1
    params:
        "-p vcf",
    log:
        "logs/{tool}/{subcommand}/tabix/index/{sample}.log",
    wrapper:
        str(wrapper_prefix / "bio" / "tabix" / "index")


""" 
This rule compress a provided VCF file with pbgzip
"""


rule pbgzip_compress:
    input:
        "{tool}/{subcommand}/{sample}.vcf",
    output:
        temp("{tool}/{subcommand}/{sample}.vcf.gz"),
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_1h_per_attempt,
        tmpdir=tmp,
    retries: 1
    params:
        "",
    log:
        "logs/{tool}/{subcommand}/pbgzip/{sample}.log",
    wrapper:
        str(wrapper_prefix / "bio" / "bcftools" / "view")
