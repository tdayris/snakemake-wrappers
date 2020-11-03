.. _`fastq_screen_indexer`:

FASTQ_SCREEN_INDEXER
====================

Build indexes and congif file for FastQ Screen



Used wrappers
---------------------


* bio/reference/ensembl-sequence

* bio/bowtie2/build




Example
-------

This meta-wrapper can be used in the following way:

.. code-block:: python

    build_release_organism = [
        "GRCh38.99.homo_sapiens",
        "GRCm38.99.mus_musculus"
    ]


    rule get_genome:
        output:
            "refs/{build}.{release}.{organism}.dna.fasta"
        params:
            species="{organism}",
            datatype="{datatype}",
            build="{build}",
            release="{release}"
        log:
            "logs/get_genome/{build}.{release}.{organism}.{datatype}.log"
        cache: True  # save space and time with between workflow caching (see docs)
        wrapper:
            "0.66.0-321-g890d65f7/bio/reference/ensembl-sequence"


Note that input, output and log file paths can be chosen freely.
When running with

.. code-block:: bash

    snakemake --use-conda

the software dependencies will be automatically deployed into an isolated environment before execution.




Authors
-------


* Thibault Dayris



Code
----


* bio/reference/ensembl-sequence

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




* bio/bowtie2/build

.. code-block:: python

    """Snakemake wrapper for bowtie2 build"""

    __author__ = "Thibault Dayris"
    __copyright__ = "Copyright 2020"
    __email__ = "koester@jimmy.harvard.edu"
    __license__ = "MIT"


    from snakemake.shell import shell
    from os.path import splitext

    extra = snakemake.params.get("extra", "")
    log = snakemake.log_fmt_shell(stdout=True, stderr=True)

    input = ""
    if "fasta" in snakemake.input.keys():
        input = "-f {}".format(snakemake.input["fasta"])
    elif "fasta" in snakemake.params.keys():
        input = "-c {}".format(snakemake.params["fasta"])
    else:
        raise ValueError(
            "Input sequence could not be found."
        )

    prefix = "bwt2_index"
    if "prefix" in snakemake.params.keys():
        prefix = snakemake.params["prefix"]


    shell(
        " bowtie2-build "
        " {input} "
        " {prefix} "
        " --threads {snakemake.threads} "
        " {extra} "
        " {log} "
    )




