wildcard_constraints:
    datatype=r"dna|cdna|cds",


rule reference_ensembl_sequence_download:
    """download genome sequence"""
    output:
        temp(
            "<temp>/reference_ensembl_sequence_download/{species}.{build}.{release}.{datatype}.fasta"
        ),
    log:
        "<logs>/reference_ensembl_sequence_download/{species}.{build}.{release}.{datatype}.log",
    benchmark:
        "<benchmarks>/curl/reference_ensembl_sequence_download_{species}.{build}.{release}.{datatype}.tsv"
    threads: 1
    params:
        species="{species}",
        datatype="{datatype}",
        build="{build}",
        release="{release}",
    wrapper:
        "v9.12.0/bio/reference/ensembl-sequence"


rule pyfaidx_filter:
    """filter out non-cannonical chromosomes"""
    input:
        fasta="<temp>/reference_ensembl_sequence_download/{species}.{build}.{release}.{datatype}.fasta",
    output:
        "<resources>/{species}.{build}.{release}/sequences/{species}.{build}.{release}.{datatype}.fasta",
    log:
        "<logs>/pyfaidx_filter/{species}.{build}.{release}.{datatype}.log",
    benchmark:
        "<benchmarks>/pyfaidx/filter_{species}.{build}.{release}.{datatype}.tsv"
    threads: 1
    params:
        extra="",
        regions=lambda wildcards: config[str(wildcards.species)][
            "cannonical_chromosomes"
        ],
    wrapper:
        "v9.4.2/bio/pyfaidx"


rule samtools_faidx:
    """index fasta sequence with samtools"""
    input:
        "<resources>/{species}.{build}.{release}/sequences/{species}.{build}.{release}.{datatype}.fasta",
    output:
        "<resources>/{species}.{build}.{release}/sequences/{species}.{build}.{release}.{datatype}.fasta.fai",
    log:
        "<logs>/samtools_faidx/{species}.{build}.{release}.{datatype}.log",
    benchmark:
        "<benchmarks>/samtools_faidx/{species}.{build}.{release}.{datatype}.tsv"
    threads: 1
    params:
        extra="",
    wrapper:
        "v9.14.0/bio/samtools/faidx"


rule picard_create_sequence_dictionary:
    """create a genome sequence dictionary with Picard"""
    input:
        "<resources>/{species}.{build}.{release}/sequences/{species}.{build}.{release}.{datatype}.fasta",
    output:
        "<resources>/{species}.{build}.{release}/sequences/{species}.{build}.{release}.{datatype}.dict",
    log:
        "<logs>/picard_create_sequence_dictionary/{species}.{build}.{release}.{datatype}.log",
    benchmark:
        "<benchmarks>/picard_create_sequence_dictionary/{species}.{build}.{release}.{datatype}.tsv"
    threads: 1
    resources:
        mem_mb=2_000,
    params:
        extra="",
    wrapper:
        "v9.4.2/bio/picard/createsequencedictionary"


rule fasta_to_twobit_convert:
    """convert genome sequence from fasta to twobit"""
    input:
        "<resources>/{species}.{build}.{release}/sequences/{species}.{build}.{release}.{datatype}.fasta",
    output:
        "<resources>/{species}.{build}.{release}/sequences/{species}.{build}.{release}.{datatype}.2bit",
    log:
        "<logs>/fasta_to_twobit_convert/{species}.{build}.{release}.{datatype}.log",
    benchmark:
        "<benchmarks>/fasta_to_twobit/convert_{species}.{build}.{release}.{datatype}.tsv"
    params:
        "",
    wrapper:
        "v7.1.0/bio/ucsc/faToTwoBit"


rule xsv_select_chrom_sizes:
    """extract chromosome sizes from fasta index"""
    input:
        table="<resources>/{species}.{build}.{release}/sequences/{species}.{build}.{release}.{datatype}.fasta.fai",
    output:
        temp("<temp>/xsv_select_chrom_sizes/{species}.{build}.{release}.{datatype}.csv"),
    log:
        "<logs>/xsv_select_chrom_sizes/{species}.{build}.{release}.{datatype}.log",
    benchmark:
        "<benchmarks>/xsv/select_chrom_sizes_{species}.{build}.{release}.{datatype}.tsv"
    threads: 1
    params:
        subcommand="select",
        extra="--no-headers --delimiter $'\\t' 1,2",
    wrapper:
        "v3.4.0/utils/xsv"


use rule xsv_select_chrom_sizes as xsv_format_chrom_sizes with:
    """format chromosome sizes from csv to tsv"""
    input:
        table="<temp>/xsv_select_chrom_sizes/{species}.{build}.{release}.{datatype}.csv",
    output:
        "<resources>/{species}.{build}.{release}/sequences/{species}.{build}.{release}.{datatype}.chrom_sizes.tsv",
    log:
        "<logs>/xsv_format_chrom_sizes/{species}.{build}.{release}.{datatype}.log",
    benchmark:
        "<benchmarks>/xsv/format_chrom_sizes_{species}.{build}.{release}.{datatype}.tsv"
    params:
        subcommand="fmt",
        extra="--out-delimiter $'\\t'",
