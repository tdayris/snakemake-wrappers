.. _`fastq_screen_indexer`:

FASTQ_SCREEN_INDEXER
====================

Build indexes and congif file for FastQ Screen


Example
-------

This meta-wrapper can be used in the following way:

.. code-block:: python

    from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
    HTTP = HTTPRemoteProvider()

    ensembl_build_release_organism = [
        ["GRCh38", "99", "homo_sapiens"],
        ["GRCm38", "99", "mus_musculus"],
        ["BDGP6.28", "99", "drosophila_melanogaster"],
        ["WBcel235", "99", "caenorhabditis_elegans"],
        ["ICSASG_v2", "99", "salmo_salar"],
        ["R64-1-1", "99", "saccharomyces_cerevisiae"],
        ["Rnor_6.0", "99", "rattus_norvegicus"]
    ]

    wildcard_constraints:
        release = r"|".join([e[1] for e in ensembl_build_release_organism]),
        organism = r"|".join([e[2] for e in ensembl_build_release_organism]),
        build = r"|".join([e[0] for e in ensembl_build_release_organism])


    rule target:
        input:
            expand(
                "index/ensembl/{fasta}.dna.1.bt2",
                fasta = [".".join(e) for e in ensembl_build_release_organism]
            )


    rule get_genome_from_ensembl:
        output:
            "refs/ensembl/{build}.{release}.{organism}.{datatype}.fasta"
        params:
            species="{organism}",
            datatype="{datatype}",
            build="{build}",
            release="{release}"
        log:
            "logs/ensembl/{build}.{release}.{organism}.{datatype}.log"
        cache: True  # save space and time with between workflow caching (see docs)
        wrapper:
            "bio/reference/ensembl-sequence"


    rule download_tair10:
        input:
            HTTP.remote(
                "https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas",
                keep_local=True
            )
        output:
            "refs/TAIR/TAIR.10.arabidopsis_thaliana.dna.fasta"
        log:
            "logs/tair/TAIR.10.arabidopsis_thaliana.chr_all.log"
        cache: True
        shell:
            "mv {input} {output} 2> {log}"


    rule download_phix_nbci:
        input:
            HTTP.remote(
                "https://ftp.ncbi.nlm.nih.gov/genomes/Viruses/enterobacteria_phage_phix174_sensu_lato_uid14015/NC_001422.fna",
                keep_local=True
            )
        output:
            "refs/ncbi/NC_001422.174.enterobacteria_phage_phix.dna.fasta"
        log:
            "logs/ncbi/NC_001422.174.enterobacteria_phage_phix.dna.log"
        cache: True
        shell:
            "mv {input} {output} 2> {log}"


    rule bowtie2_build:
        input:
            "refs/{source}/{build}.{release}.{organism}.dna.fasta"
        output:
            multiext(
                "index/{source}/{build}.{release}.{organism}.dna",
                ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2",
            )
        log:
            "logs/bwt2_build/{source}.{build}.{release}.{organism}.{datatype}.log"
        cache: True
        wrapper:
            "bio/bowtie2/build"


Note that input, output and log file paths can be chosen freely.
When running with

.. code-block:: bash

    snakemake --use-conda

the software dependencies will be automatically deployed into an isolated environment before execution.



Used wrappers
---------------------


* ``bio/reference/ensembl-sequence``

* ``bio/bowtie2/build``







Authors
-------


* Thibault Dayris



Code
----


* ``bio/reference/ensembl-sequence``

.. code-block:: python

    __author__ = "Johannes Köster"
    __copyright__ = "Copyright 2019, Johannes Köster"
    __email__ = "johannes.koester@uni-due.de"
    __license__ = "MIT"

    import subprocess as sp
    import sys
    from itertools import product
    from snakemake.shell import shell

    species = snakemake.params.species.lower()
    release = int(snakemake.params.release)
    build = snakemake.params.build

    branch = ""
    if release >= 81 and build == "GRCh37":
        # use the special grch37 branch for new releases
        branch = "grch37/"

    log = snakemake.log_fmt_shell(stdout=False, stderr=True)

    spec = ("{build}" if int(release) > 75 else "{build}.{release}").format(
        build=build, release=release
    )

    suffixes = ""
    datatype = snakemake.params.get("datatype", "")
    chromosome = snakemake.params.get("chromosome", "")
    if datatype == "dna":
        if chromosome:
            suffixes = ["dna.chromosome.{}.fa.gz".format(chromosome)]
        else:
            suffixes = ["dna.primary_assembly.fa.gz", "dna.toplevel.fa.gz"]
    elif datatype == "cdna":
        suffixes = ["cdna.all.fa.gz"]
    elif datatype == "cds":
        suffixes = ["cds.all.fa.gz"]
    elif datatype == "ncrna":
        suffixes = ["ncrna.fa.gz"]
    elif datatype == "pep":
        suffixes = ["pep.all.fa.gz"]
    else:
        raise ValueError("invalid datatype, must be one of dna, cdna, cds, ncrna, pep")

    if chromosome:
        if not datatype == "dna":
            raise ValueError(
                "invalid datatype, to select a single chromosome the datatype must be dna"
            )

    success = False
    for suffix in suffixes:
        url = "ftp://ftp.ensembl.org/pub/{branch}release-{release}/fasta/{species}/{datatype}/{species_cap}.{spec}.{suffix}".format(
            release=release,
            species=species,
            datatype=datatype,
            spec=spec.format(build=build, release=release),
            suffix=suffix,
            species_cap=species.capitalize(),
            branch=branch,
        )

        try:
            shell("curl -sSf {url} > /dev/null 2> /dev/null")
        except sp.CalledProcessError:
            continue

        shell("(curl -L {url} | gzip -d > {snakemake.output[0]}) {log}")
        success = True
        break

    if not success:
        print(
            "Unable to download requested sequence data from Ensembl. "
            "Did you check that this combination of species, build, and release is actually provided?",
            file=sys.stderr,
        )
        exit(1)




* ``bio/bowtie2/build``

.. code-block:: python

    __author__ = "Daniel Standage"
    __copyright__ = "Copyright 2020, Daniel Standage"
    __email__ = "daniel.standage@nbacc.dhs.gov"
    __license__ = "MIT"


    from snakemake.shell import shell

    extra = snakemake.params.get("extra", "")
    log = snakemake.log_fmt_shell(stdout=True, stderr=True)
    indexbase = snakemake.output[0].replace(".1.bt2", "")
    shell(
        "bowtie2-build --threads {snakemake.threads} {snakemake.params.extra} "
        "{snakemake.input.reference} {indexbase}"
    )




