wildcard_constraints:
    datatype=r"dna|cdna|cds",


rule download_genome_sequences_reference_ensembl:
    """download genome sequence"""
    output:
        temp(
            "temp/download_genome_sequences_reference_ensembl/{species}.{build}.{release}.{datatype}.fasta"
        ),
    log:
        "logs/download_genome_sequences_reference_ensembl/{species}.{build}.{release}.{datatype}.log",
    benchmark:
        "benchmarks/curl/download_genome_sequences_reference_ensembl_{species}.{build}.{release}.{datatype}.tsv"
    group:
        "download_and_filter_genome_sequence"
    threads: 1
    params:
        species="{species}",
        datatype="{datatype}",
        build="{build}",
        release="{release}",
    wrapper:
        f"{config['repo']}/bio/reference/ensembl-sequence"


rule download_genome_sequences_pyfaidx_filter:
    """filter out non-cannonical chromosomes"""
    input:
        fasta="temp/download_genome_sequences_reference_ensembl/{species}.{build}.{release}.{datatype}.fasta",
    output:
        "resources/{species}.{build}.{release}/sequences/{species}.{build}.{release}.{datatype}.fasta",
    log:
        "logs/download_genome_sequences_pyfaidx_filter/{species}.{build}.{release}.{datatype}.log",
    benchmark:
        "benchmarks/pyfaidx/download_genome_sequences_pyfaidx_filter_{species}.{build}.{release}.{datatype}.tsv"
    group:
        "download_and_filter_genome_sequence"
    threads: 1
    params:
        extra="",
        regions=lambda wildcards: config[str(wildcards.species)][
            "cannonical_chromosomes"
        ],
    wrapper:
        f"{config['repo']}/bio/pyfaidx"


rule download_genome_sequences_samtools_faidx:
    """index fasta sequence with samtools"""
    input:
        "resources/{species}.{build}.{release}/sequences/{species}.{build}.{release}.{datatype}.fasta",
    output:
        "resources/{species}.{build}.{release}/sequences/{species}.{build}.{release}.{datatype}.fasta.fai",
    log:
        "logs/download_genome_sequences_samtools_faidx/{species}.{build}.{release}.{datatype}.log",
    benchmark:
        "benchmarks/download_genome_sequences_samtools_faidx/{species}.{build}.{release}.{datatype}.tsv"
    group:
        "index_genome_sequence"
    threads: 1
    params:
        extra="",
    wrapper:
        f"{config['repo']}/bio/samtools/faidx"


rule download_genome_sequences_picard_create_sequence_dictionary:
    """create a genome sequence dictionary with Picard"""
    input:
        "resources/{species}.{build}.{release}/sequences/{species}.{build}.{release}.{datatype}.fasta",
    output:
        "resources/{species}.{build}.{release}/sequences/{species}.{build}.{release}.{datatype}.dict",
    log:
        "logs/download_genome_sequences_picard_create_sequence_dictionary/{species}.{build}.{release}.{datatype}.log",
    benchmark:
        "benchmarks/download_genome_sequences_picard_create_sequence_dictionary/{species}.{build}.{release}.{datatype}.tsv"
    threads: 1
    resources:
        mem_mb=2_000,
    params:
        extra="",
    wrapper:
        f"{config['repo']}/bio/picard/createsequencedictionary"


rule download_genome_sequences_fasta_to_twobit_convert:
    """convert genome sequence from fasta to twobit"""
    input:
        "resources/{species}.{build}.{release}/sequences/{species}.{build}.{release}.{datatype}.fasta",
    output:
        "resources/{species}.{build}.{release}/sequences/{species}.{build}.{release}.{datatype}.2bit",
    log:
        "logs/download_genome_sequences_fasta_to_twobit_convert/{species}.{build}.{release}.{datatype}.log",
    benchmark:
        "benchmarks/fasta_to_twobit/download_genome_sequences_fasta_to_twobit_convert_{species}.{build}.{release}.{datatype}.tsv"
    params:
        "",
    wrapper:
        f"{config['repo']}/bio/ucsc/faToTwoBit"


rule download_genome_sequences_xan_extract_chrom_sizes:
    """extract chromosome sizes from fasta index"""
    input:
        data="resources/{species}.{build}.{release}/sequences/{species}.{build}.{release}.{datatype}.fasta.fai",
    output:
        "resources/{species}.{build}.{release}/sequences/{species}.{build}.{release}.{datatype}.chrom_sizes.tsv",
    log:
        "logs/download_genome_sequences_xan_extract_chrom_sizes/{species}.{build}.{release}.{datatype}.log",
    benchmark:
        "benchmarks/xan/download_genome_sequences_xan_extract_chrom_sizes_{species}.{build}.{release}.{datatype}.tsv"
    group:
        "index_genome_sequence"
    threads: 1
    params:
        expression="select 1,2 | fmt --tabs",
        extra="--tee",
    wrapper:
        f"{config['repo']}/utils/xan/run"


rule download_genome_sequences_sed_replace_genomes_data_in_report:
    """replace snakemake patterns for generic report"""
    input:
        workflow.source_path(
            "resources/reports/material_and_methods_genome_sequence.rst"
        ),
    output:
        temp(
            "temp/download_genome_sequences_sed_replace_genomes_data_in_report/material_and_methods_{species}_{build}_{release}.rst"
        ),
    log:
        "logs/download_genome_sequences_sed_replace_genomes_data_in_report/material_and_methods_{species}_{build}_{release}.log",
    benchmark:
        "benchmarks/sed/download_genome_sequences_sed_replace_genomes_data_in_report_{species}_{build}_{release}.tsv"
    group:
        "genome_sequence_material_and_methods"
    threads: 1
    params:
        extra="",
        expr="s| params.species |{species}|g;s| params.build |{build}|g;s| params.release |{release}|g;s|\\{\\{||g;s|\\}\\}||g",
    wrapper:
        f"{config['repo']}/utils/sed"


rule download_genome_sequences_material_and_methods:
    """provide information on genome sequence preparation methods"""
    input:
        rst="temp/download_genome_sequences_sed_replace_genomes_data_in_report/material_and_methods_{species}_{build}_{release}.rst",
    output:
        html="resources/{species}.{build}.{release}/sequences/Genome_sequences_material_and_methods.html",
    log:
        "logs/download_genome_sequences_material_and_methods/{species}.{build}.{release}.log",
    benchmark:
        "benchmarks/rst2html/download_genome_sequences_material_and_methods_{species}.{build}.{release}.tsv"
    group:
        "genome_sequence_material_and_methods"
    threads: 1
    params:
        extra="--title='Material and methods for genome sequence preparation' --embed-stylesheet",
    wrapper:
        f"{config['repo']}/utils/docutils/rst2html"
