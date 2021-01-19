.. _`Download Reference`:

DOWNLOAD REFERENCE
==================

Download sequence and annotations from Ensembl


Example
-------

This meta-wrapper can be used by integrating the following into your workflow:

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
            "v0.69.0-406-g1e9e192ca/bio/reference/ensembl-sequence"


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
            "v0.69.0-406-g1e9e192ca/bio/reference/ensembl-annotation"


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
            "v0.69.0-406-g1e9e192ca/bio/samtools/faidx"


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
            "v0.69.0-406-g1e9e192ca/bio/picard/createsequencedictionary"


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
            "v0.69.0-406-g1e9e192ca/bio/reference/ensembl-variation"

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

* :ref:`bio/reference/ensembl-annotation`

* :ref:`bio/reference/ensembl-variation`

* :ref:`bio/samtools/faidx`

* :ref:`bio/picard/createsequencedictionary`


Please refer to each wrapper in above list for additional configuration parameters and information about the executed code.






Notes
-----

Do not forget to cache these downloads!

The samtools index step is here to include genome intervals in the VCF index. By doing so, the VCF is compatible with GATK for variant calling.




Authors
-------


* Thibault Dayris

