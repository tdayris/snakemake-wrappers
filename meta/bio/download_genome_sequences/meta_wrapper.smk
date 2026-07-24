rule reference_ensembl_sequence_download:
    """download genome sequence"""
    output:
        temp(
            "<tmp>/reference_ensembl_sequence_download/{species}.{build}.{release}.{datatype}.fasta"
        ),
    log:
        "<log>/reference_ensembl_sequence_download/{species}.{build}.{release}.{datatype}.log",
    benchmark:
        "<benchmark>/reference_ensembl_sequence_download/{species}.{build}.{release}.{datatype}.tsv"
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
        fasta="<tmp>/reference_ensembl_sequence_download/{species}.{build}.{release}.{datatype}.fasta",
    output:
        "<reference>/{species}.{build}.{release}/sequences/{species}.{build}.{release}.{datatype}.fasta",
    log:
        "<log>/pyfaidx_filter/{species}.{build}.{release}.{datatype}.log",
    benchmark:
        "<benchmark>/pyfaidx_filter/{species}.{build}.{release}.{datatype}.tsv"
    threads: 1
    params:
        extra="",
        regions=config["contigs"],
    wrapper:
        "v9.4.2/bio/pyfaidx"


rule samtools_faidx:
    """index fasta sequence with samtools"""
    input:
        "<reference>/{species}.{build}.{release}/sequences/{species}.{build}.{release}.{datatype}.fasta",
    output:
        "<reference>/{species}.{build}.{release}/sequences/{species}.{build}.{release}.{datatype}.fasta.fai",
    log:
        "<log>/samtools_faidx/{species}.{build}.{release}.{datatype}.log",
    benchmark:
        "<benchmark>/samtools_faidx/{species}.{build}.{release}.{datatype}.tsv"
    threads: 1
    params:
        extra="",
    wrapper:
        "v9.14.0/bio/samtools/faidx"


rule picard_create_sequence_dictionary:
    """create a genome sequence dictionary with Picard"""
    input:
        "<reference>/{species}.{build}.{release}/sequences/{species}.{build}.{release}.{datatype}.fasta",
    output:
        "<reference>/{species}.{build}.{release}/sequences/{species}.{build}.{release}.{datatype}.dict",
    log:
        "<log>/picard_create_sequence_dictionary/{species}.{build}.{release}.{datatype}.log",
    benchmark:
        "<benchmark>/picard_create_sequence_dictionary/{species}.{build}.{release}.{datatype}.tsv"
    threads: 1
    params:
        extra="",
    wrapper:
        "v9.4.2/bio/picard/createsequencedictionary"


rule fasta_to_twobit_convert:
    """convert genome sequence from fasta to twobit"""
    input:
        "<reference>/{species}.{build}.{release}/sequences/{species}.{build}.{release}.{datatype}.fasta",
    output:
        "<reference>/{species}.{build}.{release}/sequences/{species}.{build}.{release}.{datatype}.2bit",
    log:
        "<log>/fasta_to_twobit_convert/{species}.{build}.{release}.{datatype}.log",
    benchmark:
        "<benchmark>/fasta_to_twobit_convert/{species}.{build}.{release}.{datatype}.tsv",
    params:
        "",
    wrapper:
        "v7.1.0/bio/ucsc/faToTwoBit"


rule xsv_select_chrom_sizes:
    """extract chromosome sizes from fasta index"""
    input:
        table="<reference>/{species}.{build}.{release}/sequences/{species}.{build}.{release}.{datatype}.fasta.fai",
    output:
        temp("<tmp>/xsv_select_chrom_sizes/{species}.{build}.{release}.{datatype}.csv"),
    log:
        "<log>/xsv_select_chrom_sizes/{species}.{build}.{release}.{datatype}.log",
    benchmark:
        "<benchmark/xsv_select_chrom_sizes/{species}.{build}.{release}.{datatype}.tsv"
    threads: 1
    params:
        subcommand="select",
        extra="--no-headers --delimiter $'\\t' 1,2",
    wrapper:
        "v3.4.0/utils/xsv"


rule xsv_format_chrom_sizes:
    """format chromosome sizes from csv to tsv"""
    input:
        table="<tmp>/xsv_select_chrom_sizes/{species}.{build}.{release}.{datatype}.csv",
    output:
        "<reference>/{species}.{build}.{release}/sequences/{species}.{build}.{release}.{datatype}.chrom_sizes.tsv",
    log:
        "<log>/xsv_format_chrom_sizes/{species}.{build}.{release}.{datatype}.log"
    benchmark:
        "<benchmark>/xsv_format_chrom_sizes/{species}.{build}.{release}.{datatype}.tsv"
    params:
        subcommand="fmt",
        extra="--out-delimiter $'\\t'",
    wrapper:
        "v3.4.0/utils/xsv"
