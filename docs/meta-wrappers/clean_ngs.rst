.. _`clean_ngs`:

CLEAN_NGS
=========

Clean and control short reads with FastqScreen, Fastp and MultiQC


Example
-------

This meta-wrapper can be used by integrating the following into your workflow:

.. code-block:: python

    from snakemake.utils import min_version
    min_version("6.0")

    module fastq_screen_indexer:
        meta_wrapper: "0.71.1-451-gb2e59cf65/meta/bio/fastq_screen_indexer"


    build_release_organism_list = [
        ["ensembl", "GRCh38", "99", "homo_sapiens"],
        ["ensembl", "GRCm38", "99", "mus_musculus"],
        ["ensembl", "BDGP6.28", "99", "drosophila_melanogaster"],
        ["ensembl", "WBcel235", "99", "caenorhabditis_elegans"],
        ["ensembl", "ICSASG_v2", "99", "salmo_salar"],
        ["ensembl", "R64-1-1", "99", "saccharomyces_cerevisiae"],
        ["ensembl", "Rnor_6.0", "99", "rattus_norvegicus"],
        ["TAIR", "TAIR", "10", "arabidopsis_thaliana"],
        ["ncbi", "NC_001422", "174", "enterobacteria_phage_phix"]
    ]


    fastq_screen_config = {
        "database": {
            organism: {f"index/{source}/{build}.{release}.{organism}.dna"}
            for source, build, release, organism in build_release_organism_list
        }
    }

    use rule * from fastq_screen_indexer

    rule multiqc_report:
        input:
            fastq_screen=expand(
                "fastq_screen/{sample}.fastq_screen.{ext}",
                ext=["txt", "png"],
                sample=config.get("sample_list", ["a"])
            ),
            fastp=expand(
                "fastp/{format}/se/{sample}.fastp.{format}",
                format=["html", "json"],
                sample=config.get("sample_list", ["a"])
            )
        output:
            "multiqc/report.html"
        message: "Gathering samples quality with MultiQC"
        threads: 1
        resources:
            mem_mb=lambda wildcard, attempt: min(attempt * 2048, 10240),
            time_min=lambda wildcard, attempt: attempt * 30
        params:
            config.get("multiqc_extra", "")
        log:
            "logs/multiqc.log"
        wrapper:
            "0.71.1-453-g032eb4537/bio/multiqc"


    rule fastq_screen_analysis:
        input:
            "reads/{sample}.fastq",
            bt2_index=expand(
                "index/{source}/{reference}.dna.bt2",
                source=[i[0] for i in build_release_organism_list],
                reference=[".".join(i[1:]) for i in build_release_organism_list]
            )
        output:
            txt="fastq_screen/{sample}.fastq_screen.txt",
            png="fastq_screen/{sample}.fastq_screen.png"
        message: "Assessing quality with FastqScreen on {wildcard.sample}"
        threads: config.get("threads", 20)
        resources:
            mem_mb=lambda wildcard, attempt: min(attempt * 8192, 15360),
            time_min=lambda wildcard, attempt: attempt * 45
        params:
            fastq_screen_config=fastq_screen_config,
            subset=config.get("fastq_screen_subset", 100000),
            aligner="bowtie2",
            extra=config.get("fastq_screen_extra", "")
        log:
            "logs/fastq_screen/analysis/{sample}.log"
        wrapper:
            "0.71.1-453-g032eb4537/bio/fastq_screen"


    rule fastp_clean:
        input:
            sample=expand(
                "reads/{sample}.{stream}.fastq",
                stream=["1", "2"],
                allow_missing=True
            ),
        output:
            trimmed=expand(
                "fastp/trimmed/pe/{sample}.{stream}.fastq",
                "fastp/trimmed/pe/{sample}.{stream}.fastq",
                stream=["1", "2"],
                allow_missing=True
            ),
            html="fastp/html/pe/{sample}.fastp.html",
            json="fastp/json/pe/{sample}.fastp.json"
        message: "Cleaning {wildcard.sample} with Fastp"
        threads: 1
        resources:
            mem_mb=lambda wildcard, attempt: min(attempt * 4096, 15360),
            time_min=lambda wildcard, attempt: attempt * 45
        params:
            adapters=config.get("fastp_adapters", ""),
            extra=config.get("fastp_extra", "")
        log:
            "logs/fastp/{sample}.log"
        wrapper:
            "0.71.1-453-g032eb4537/bio/fastp"

Note that input, output and log file paths can be chosen freely, as long as the dependencies between the rules remain as listed here.
For additional parameters in each individual wrapper, please refer to their corresponding documentation (see links below).

When running with

.. code-block:: bash

    snakemake --use-conda

the software dependencies will be automatically deployed into an isolated environment before execution.



Used wrappers
---------------------

The following individual wrappers are used in this meta-wrapper:


* :ref:`bio/fastp`

* :ref:`bio/fastq_screen`


Please refer to each wrapper in above list for additional configuration parameters and information about the executed code.






Notes
-----

MultiQC sums up all the results.




Authors
-------


* Thibault Dayris

