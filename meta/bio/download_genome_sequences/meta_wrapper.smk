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
        f"{config['repo']}/bio/reference/ensembl-sequence"


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
        f"{config['repo']}/bio/pyfaidx"


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
        f"{config['repo']}/bio/samtools/faidx"


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
        f"{config['repo']}/bio/picard/createsequencedictionary"


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
        f"{config['repo']}/bio/ucsc/faToTwoBit"


rule xan_extract_chrom_sizes:
    """extract chromosome sizes from fasta index"""
    input:
        table="<resources>/{species}.{build}.{release}/sequences/{species}.{build}.{release}.{datatype}.fasta.fai",
    output:
        "<resources>/{species}.{build}.{release}/sequences/{species}.{build}.{release}.{datatype}.chrom_sizes.tsv",
    log:
        "<logs>/xan_extract_chrom_sizes/{species}.{build}.{release}.{datatype}.log",
    benchmark:
        "<benchmarks>/xan/select_chrom_sizes_{species}.{build}.{release}.{datatype}.tsv"
    threads: 1
    params:
        expression="select 1,2 | fmt --tabs",
        extra="--tee",
    wrapper:
        f"{config['repo']}/utils/xan/run"


rule material_and_methods_genome_sequences:
    """provide information on genome sequence preparation methods"""
    input:
        rst=workflow.source_path("resources/reports/material_and_methods.rst"),
    output:
        html="<resources>/{species}.{build}.{release}/sequence/Genome_sequences_material_and_methods.html",
    log:
        "<logs>/material_and_methods_genome_sequences/{species}.{build}.{release}.log",
    benchmark:
        "<benchmarks>/rst2html/material_and_methods_genome_sequences_{species}.{build}.{release}.tsv"
    threads: 1
    params:
        extra="--title='Material and methods for genome sequence preparation' --embed-stylesheet",
    wrapper:
        f"{config['repo']}/utils/docutils/rst2html"
