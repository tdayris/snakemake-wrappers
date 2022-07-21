
module compress_index_vcf_meta:
    snakefile:
        str(
            worflow_source_dir
            / ".."
            / ".."
            / "meta"
            / "bio"
            / "compress_index_vcf"
            / "test"
            / "Snakefile"
        )
    config:
        config


use rule * from compress_index_vcf_meta


use rule tabix_index from compress_index_vcf_meta as snp_indel_tabix_index with:
    input:
        "{tool}/{subcommand}/{sample}.{content}.vcf.gz",
    output:
        "{tool}/{subcommand}/{sample}.{content}.vcf.gz.tbi",
    message:
        "Indexing {wildcards.sample} (from somatic varscan "
        "{wildcards.content}) with tabix."
    log:
        "logs/{tool}/{subcommand}/tabix/index/{sample}.{content}.log",


use rule pbgzip_compress from compress_index_vcf_meta as si_pbgzip with:
    input:
        "{tool}/{subcommand}/{sample}.{content}.vcf",
    output:
        "{tool}/{subcommand}/{sample}.{content}.vcf.gz",
    message:
        "Compressnig {wildcards.sample} (from somatic varscan "
        "{wildcards.content}) with pbgzip."
    log:
        "logs/{tool}/{subcommand}/pgbzip/varcsanc2/{sample}.{content}.log",
