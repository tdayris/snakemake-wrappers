.. _`fastq_screen_indexer`:

FASTQ_SCREEN_INDEXER
====================

Build indexes and congif file for FastQ Screen


Example
-------

This meta-wrapper can be used by integrating the following into your workflow:

.. code-block:: python

    from snakemake.utils import min_version
    min_version("6.0")

    from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
    HTTP = HTTPRemoteProvider()

    localrules: get_genome, download_tair10, download_phix_nbci


    build_release_organism = [
        ["GRCh38", "99", "homo_sapiens", "ensembl"],
        ["GRCm38", "99", "mus_musculus", "ensembl"],
        ["BDGP6.28", "99", "drosophila_melanogaster", "ensembl"],
        ["WBcel235", "99", "caenorhabditis_elegans", "ensembl"],
        ["ICSASG_v2", "99", "salmo_salar", "ensembl"],
        ["R64-1-1", "99", "saccharomyces_cerevisiae", "ensembl"],
        ["Rnor_6.0", "99", "rattus_norvegicus", "ensembl"],
        ["TAIR", "10", "arabidopsis_thaliana", "TAIR"],
        ["NC_001422", "174", "enterobacteria_phage_phix", "ncbi"]
    ]

    wildcard_constraints:
        release = r"|".join([e[1] for e in build_release_organism]),
        organism = r"|".join(set([e[2] for e in build_release_organism])),
        build = r"|".join([e[0] for e in build_release_organism]),
        source = r"|".join(set([e[3] for e in build_release_organism]))


    rule target:
        input:
            [
                f"refs/{source}/{build}.{release}.{organism}.dna.fasta"
                #f"index/{source}/{build}.{release}.{organism}.dna.1.bt2"
                for build, release, organism, source in build_release_organism
            ]


    rule get_genome:
        output:
            "refs/ensembl/{build}.{release}.{organism}.{datatype}.fasta"
        threads: 1
        resources:
            mem_mb=lambda wildcard, attempt: min(attempt * 512, 2048),
            time_min=lambda wildcard, attempt: attempt * 120
        params:
            species="{organism}",
            datatype="{datatype}",
            build="{build}",
            release="{release}"
        log:
            "logs/get_genome/{build}.{release}.{organism}.{datatype}.log"
        cache: True  # save space and time with between workflow caching (see docs)
        wrapper:
            "0.75.0-757-g5432694f9/bio/reference/ensembl-sequence"


    rule download_tair10:
        input:
            HTTP.remote(
                "https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas",
                keep_local=True
            )
        output:
            "refs/TAIR/TAIR.10.arabidopsis_thaliana.dna.fasta"
        cache: True
        threads: 1
        resources:
            mem_mb=lambda wildcard, attempt: min(attempt * 512, 2048),
            time_min=lambda wildcard, attempt: attempt * 120
        log:
            "logs/tair/TAIR.10.arabidopsis_thaliana.chr_all.log"
        shell:
            "mv {input} {output} 2> {log}"


    use rule download_tair10 as download_phix_nbci with:
        input:
            HTTP.remote(
                "https://ftp.ncbi.nlm.nih.gov/genomes/Viruses/enterobacteria_phage_phix174_sensu_lato_uid14015/NC_001422.fna",
                keep_local=True
            )
        output:
            "refs/ncbi/NC_001422.174.enterobacteria_phage_phix.dna.fasta"
        log:
            "logs/ncbi/NC_001422.174.enterobacteria_phage_phix.dna.log"


    rule bowtie2_build:
        input:
            reference = "refs/{source}/{build}.{release}.{organism}.dna.fasta"
        output:
            multiext(
                "index/{source}/{build}.{release}.{organism}.dna",
                ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2",
            )
        cache: True
        threads: config.get("threads", 10)
        resources:
            mem_mb=lambda wildcard, attempt: min(attempt * 4096, 8192),
            time_min=lambda wildcard, attempt: attempt * 45
        log:
            "logs/bwt2_build/{source}.{build}.{release}.{organism}.dna.log"
        params:
            extra=""
        wrapper:
            "0.75.0-757-g5432694f9/bio/bowtie2/build"

Note that input, output and log file paths can be chosen freely, as long as the dependencies between the rules remain as listed here.
For additional parameters in each individual wrapper, please refer to their corresponding documentation (see links below).

When running with

.. code-block:: bash

    snakemake --use-conda

the software dependencies will be automatically deployed into an isolated environment before execution.



Used wrappers
---------------------

The following individual wrappers are used in this meta-wrapper:


* :ref:`bio/reference/ensembl-sequence`

* :ref:`bio/bowtie2/build`


Please refer to each wrapper in above list for additional configuration parameters and information about the executed code.







Authors
-------


* Thibault Dayris

