"""
001::samtools_index_genome Independent, optional rule

This rule is executed if and only if no .fasta.fai file
is provided in configfile.
"""
rule samtools_index_genome:
    input:
        config[genome_id]["fasta"]
    output:
        fai_file
    threads: 1
    resources:
        mem_mb=get_768mb_per_gb,
        time_min=get_10_minutes_per_gb,
        tmpdir="tmp"
    params:
        extra="",
    log:
        "logs/samtools/faidx/genome.log"
    wrapper:
        "bio/samtools/faidx"


"""
001::picard_createsequencedictionary Independent, optional rule

This rule makes the fasta sequence dictionary
in case user did not provided any in the configfile
"""
rule picard_createsequencedictionary:
    input:
        config[genome_id]["fasta"]
    output:
        dict_file
    threads: 1
    resources:
        mem_mb=get_768mb_per_gb,
        time_min=get_10_minutes_per_gb,
        tmpdir="tmp"
    params:
        extra=(
            f"--GENOME_ASSEMBLY {config[genome_id]['genome_assembly']} "
            f"--SPECIES {config[genome_id]['species']}"
        ),
    log:
        "logs/picard/createsequencedictionary/genome.log"
    wrapper:
        "bio/picard/createsequencedictionary"


"""
001::tabix_index Independent, optional rule

In case of missing known variants index, this
rule tab indexes the input VCF file
"""
rule tabix_index_known_variants:
    input:
        config[genome_id]["vcf"]
    output:
        tbi_file
    threads: 1
    resources:
        mem_mb=get_1p5gb_per_gb,
        time_min=get_10_minutes_per_gb,
        tmpdir="tmp"
    log:
        "logs/tabix/index/known_variants.log"
    params:
        "-p vcf"
    wrapper:
        "bio/tabix"