.. _`Download Reference`:

DOWNLOAD REFERENCE
==================

Download sequence and annotations from Ensembl



Used wrappers
---------------------


* bio/reference/ensembl-sequence

* bio/reference/ensembl-annotation

* bio/reference/ensembl-variation

* bio/samtools/faidx

* bio/picard/createsequencedictionary




Example
-------

This meta-wrapper can be used in the following way:

.. code-block:: python

    fasta_datatype = ["dna", "cdna", "ncrna"]
    build_release_organism = [
        "GRCh38.99.homo_sapiens",
        "GRCm38.99.mus_musculus"
    ]

    rule all:
        input:
            fasta = expand(
                "refs/{build_release_organism}.{datatype}.fasta",
                build_release_organism=build_release_organism,
                datatype=fasta_datatype
            ),
            fasta_index = expand(
                "refs/{build_release_organism}.{datatype}.fasta.fai",
                build_release_organism=build_release_organism,
                datatype=fasta_datatype
            ),
            fasta_dict = expand(
                "refs/{build_release_organism}.{datatype}.dict",
                build_release_organism=build_release_organism,
                datatype=fasta_datatype
            ),
            gtf = expand(
                "refs/{build_release_organism}.gtf",
                build_release_organism=build_release_organism
            ),
            vcf = expand(
                "refs/{build_release_organism}.all.vcf.gz",
                build_release_organism=build_release_organism
            )


    rule get_genome:
        output:
            "refs/{build}.{release}.{organism}.{datatype}.fasta"
        params:
            species="{organism}",
            datatype="{datatype}",
            build="{build}",
            release="{release}"
        log:
            "logs/get_genome/{build}.{release}.{organism}.{datatype}.log"
        cache: True  # save space and time with between workflow caching (see docs)
        wrapper:
            "0.66.0-280-g9f0281de/bio/reference/ensembl-sequence"


    rule get_annotation:
        output:
            "refs/{build}.{release}.{organism}.gtf"
        params:
            species="{organism}",
            release="{release}",
            build="{build}",
            fmt="gtf",
            flavor="" # optional, e.g. chr_patch_hapl_scaff, see Ensembl FTP.
        log:
            "logs/get_annotation/{build}.{release}.{organism}.log"
        cache: True  # save space and time with between workflow caching (see docs)
        wrapper:
            "0.66.0-280-g9f0281de/bio/reference/ensembl-annotation"


    rule samtools_faidx_reference:
        input:
            "refs/{build}.{release}.{organism}.{datatype}.fasta"
        output:
            "refs/{build}.{release}.{organism}.{datatype}.fasta.fai"
        params:
            "" # optional params string
        cache: True
        group: "index_fasta"
        wrapper:
            "0.66.0-280-g9f0281de/bio/samtools/faidx"


    rule create_dict:
        input:
            "refs/{build}.{release}.{organism}.{datatype}.fasta"
        output:
            "refs/{build}.{release}.{organism}.{datatype}.dict"
        log:
            "logs/picard/create_dict/{build}.{release}.{organism}.{datatype}.log"
        params:
            extra=""  # optional: extra arguments for picard.
        cache: True
        group: "index_fasta"
        wrapper:
            "0.66.0-280-g9f0281de/bio/picard/createsequencedictionary"


    rule get_variation_with_contig_lengths:
        input:
            fai="refs/{build}.{release}.{organism}.dna.fasta.fai"
        output:
            vcf="refs/{build}.{release}.{organism}.all.vcf.gz"
        params:
            species="{organism}",
            release="{release}",
            build="{build}",
            type="all" # one of "all", "somatic", "structural_variation"
        log:
            "logs/get_variation/{build}.{release}.{organism}.log"
        wrapper:
            "0.66.0-280-g9f0281de/bio/reference/ensembl-variation"


Note that input, output and log file paths can be chosen freely.
When running with

.. code-block:: bash

    snakemake --use-conda

the software dependencies will be automatically deployed into an isolated environment before execution.



Notes
-----

Do not forget to cache these downloads!

The samtools index step is here to include genome intervals in the VCF index. By doing so, the VCF is compatible with GATK for variant calling.




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
    if datatype == "dna":
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




* bio/reference/ensembl-annotation

.. code-block:: python

    __author__ = "Johannes Köster"
    __copyright__ = "Copyright 2019, Johannes Köster"
    __email__ = "johannes.koester@uni-due.de"
    __license__ = "MIT"

    import subprocess
    import sys
    from snakemake.shell import shell

    species = snakemake.params.species.lower()
    release = int(snakemake.params.release)
    fmt = snakemake.params.fmt
    build = snakemake.params.build
    flavor = snakemake.params.get("flavor", "")

    branch = ""
    if release >= 81 and build == "GRCh37":
        # use the special grch37 branch for new releases
        branch = "grch37/"

    if flavor:
        flavor += "."

    log = snakemake.log_fmt_shell(stdout=False, stderr=True)

    suffix = ""
    if fmt == "gtf":
        suffix = "gtf.gz"
    elif fmt == "gff3":
        suffix = "gff3.gz"

    url = "ftp://ftp.ensembl.org/pub/{branch}release-{release}/{fmt}/{species}/{species_cap}.{build}.{release}.{flavor}{suffix}".format(
        release=release,
        build=build,
        species=species,
        fmt=fmt,
        species_cap=species.capitalize(),
        suffix=suffix,
        flavor=flavor,
        branch=branch,
    )

    try:
        shell("(curl -L {url} | gzip -d > {snakemake.output[0]}) {log}")
    except subprocess.CalledProcessError as e:
        if snakemake.log:
            sys.stderr = open(snakemake.log[0], "a")
        print(
            "Unable to download annotation data from Ensembl. "
            "Did you check that this combination of species, build, and release is actually provided?",
            file=sys.stderr,
        )
        exit(1)




* bio/reference/ensembl-variation

.. code-block:: python

    __author__ = "Johannes Köster"
    __copyright__ = "Copyright 2019, Johannes Köster"
    __email__ = "johannes.koester@uni-due.de"
    __license__ = "MIT"

    import tempfile
    import subprocess
    import sys
    import os
    from snakemake.shell import shell
    from snakemake.exceptions import WorkflowError

    species = snakemake.params.species.lower()
    release = int(snakemake.params.release)
    build = snakemake.params.build
    type = snakemake.params.type

    if release < 98:
        print("Ensembl releases <98 are unsupported.", file=open(snakemake.log[0], "w"))
        exit(1)

    branch = ""
    if release >= 81 and build == "GRCh37":
        # use the special grch37 branch for new releases
        branch = "grch37/"

    log = snakemake.log_fmt_shell(stdout=False, stderr=True)

    if type == "all":
        if species == "homo_sapiens" and release >= 93:
            suffixes = [
                "-chr{}".format(chrom) for chrom in list(range(1, 23)) + ["X", "Y", "MT"]
            ]
        else:
            suffixes = [""]
    elif type == "somatic":
        suffixes = ["_somatic"]
    elif type == "structural_variations":
        suffixes = ["_structural_variations"]
    else:
        raise ValueError(
            "Unsupported type {} (only all, somatic, structural_variations are allowed)".format(
                type
            )
        )

    species_filename = species if release >= 91 else species.capitalize()

    urls = [
        "ftp://ftp.ensembl.org/pub/{branch}release-{release}/variation/vcf/{species}/{species_filename}{suffix}.{ext}".format(
            release=release,
            species=species,
            suffix=suffix,
            species_filename=species_filename,
            branch=branch,
            ext=ext,
        )
        for suffix in suffixes
        for ext in ["vcf.gz", "vcf.gz.csi"]
    ]
    names = [os.path.basename(url) for url in urls if url.endswith(".gz")]

    try:
        gather = "curl {urls}".format(urls=" ".join(map("-O {}".format, urls)))
        workdir = os.getcwd()
        with tempfile.TemporaryDirectory() as tmpdir:
            if snakemake.input.get("fai"):
                shell(
                    "(cd {tmpdir}; {gather} && "
                    "bcftools concat -Oz --naive {names} > concat.vcf.gz && "
                    "bcftools reheader --fai {workdir}/{snakemake.input.fai} concat.vcf.gz "
                    "> {workdir}/{snakemake.output}) {log}"
                )
            else:
                shell(
                    "(cd {tmpdir}; {gather} && "
                    "bcftools concat -Oz --naive {names} "
                    "> {workdir}/{snakemake.output}) {log}"
                )
    except subprocess.CalledProcessError as e:
        if snakemake.log:
            sys.stderr = open(snakemake.log[0], "a")
        print(
            "Unable to download variation data from Ensembl. "
            "Did you check that this combination of species, build, and release is actually provided? ",
            file=sys.stderr,
        )
        exit(1)




* bio/samtools/faidx

.. code-block:: python

    __author__ = "Michael Chambers"
    __copyright__ = "Copyright 2019, Michael Chambers"
    __email__ = "greenkidneybean@gmail.com"
    __license__ = "MIT"


    from snakemake.shell import shell


    shell("samtools faidx {snakemake.params} {snakemake.input[0]} > {snakemake.output[0]}")




* bio/picard/createsequencedictionary

.. code-block:: python

    __author__ = "Johannes Köster"
    __copyright__ = "Copyright 2018, Johannes Köster"
    __email__ = "johannes.koester@protonmail.com"
    __license__ = "MIT"


    from snakemake.shell import shell


    extra = snakemake.params.get("extra", "")
    log = snakemake.log_fmt_shell(stdout=False, stderr=True)


    memory = ""
    if "mem_mb" in snakemake.resources.keys():
        memory = "-Xmx{}M".format(snakemake.resources["mem_mb"])

    shell(
        "picard "
        "CreateSequenceDictionary "
        "{memory} "
        "{extra} "
        "R={snakemake.input[0]} "
        "O={snakemake.output[0]} "
        "{log}"
    )




