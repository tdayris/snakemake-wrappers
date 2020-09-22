.. _`star-arriba`:

STAR-ARRIBA
===========

A subworkflow for fusion detection from RNA-seq data with ``arriba``. The fusion calling is based on splice-aware, chimeric alignments done with ``STAR``. ``STAR`` is used with specific parameters to ensure optimal functionality of the ``arriba`` fusion detection, for details, see the `documentation <https://arriba.readthedocs.io/en/latest/workflow/>`_.




Used wrappers
---------------------


* bio/star/index

* bio/star/align

* bio/arriba




Example
-------

This meta-wrapper can be used in the following way:

.. code-block:: python

    rule star_index:
        input:
            fasta="resources/genome.fasta",
            annotation="resources/genome.gtf"
        output:
            directory("resources/star_genome")
        threads: 4
        params:
            extra="--sjdbGTFfile resources/genome.gtf --sjdbOverhang 100"
        log:
            "logs/star_index_genome.log"
        cache: True
        wrapper:
            "0.66.0-240-gd9cffe8c/bio/star/index"

    rule star_align:
        input:
            # use a list for multiple fastq files for one sample
            # usually technical replicates across lanes/flowcells
            fq1="reads/{sample}_R1.1.fastq",
            fq2="reads/{sample}_R2.1.fastq", #optional
            index="resources/star_genome"
        output:
            # see STAR manual for additional output files
            "star/{sample}/Aligned.out.bam",
            "star/{sample}/ReadsPerGene.out.tab"
        log:
            "logs/star/{sample}.log"
        params:
            # path to STAR reference genome index
            index="resources/star_genome",
            # specific parameters to work well with arriba
            extra="--quantMode GeneCounts --sjdbGTFfile resources/genome.gtf"
                " --outSAMtype BAM Unsorted --chimSegmentMin 10 --chimOutType WithinBAM SoftClip"
                " --chimJunctionOverhangMin 10 --chimScoreMin 1 --chimScoreDropMax 30 --chimScoreJunctionNonGTAG 0"
                " --chimScoreSeparation 1 --alignSJstitchMismatchNmax 5 -1 5 5 --chimSegmentReadGapMax 3"
        threads: 12
        wrapper:
            "0.66.0-240-gd9cffe8c/bio/star/align"

    rule arriba:
        input:
            bam="star/{sample}/Aligned.out.bam",
            genome="resources/genome.fasta",
            annotation="resources/genome.gtf"
        output:
            fusions="results/arriba/{sample}.fusions.tsv",
            discarded="results/arriba/{sample}.fusions.discarded.tsv"
        params:
            # A tsv containing identified artifacts, such as read-through fusions of neighbouring genes, see https://arriba.readthedocs.io/en/latest/input-files/#blacklist
            blacklist="arriba_blacklist.tsv",
            extra="-T -P -i 1,2" # -i describes the wanted contigs, remove if you want to use all hg38 chromosomes
        log:
            "logs/arriba/{sample}.log"
        threads: 1
        wrapper:
            "0.66.0-240-gd9cffe8c/bio/arriba"


Note that input, output and log file paths can be chosen freely.
When running with

.. code-block:: bash

    snakemake --use-conda

the software dependencies will be automatically deployed into an isolated environment before execution.




Authors
-------


* Jan Forster



Code
----


* bio/star/index

.. code-block:: python

    """Snakemake wrapper for STAR index"""

    __author__ = "Thibault Dayris"
    __copyright__ = "Copyright 2019, Dayris Thibault"
    __email__ = "thibault.dayris@gustaveroussy.fr"
    __license__ = "MIT"

    from snakemake.shell import shell
    from snakemake.utils import makedirs

    log = snakemake.log_fmt_shell(stdout=True, stderr=True)

    extra = snakemake.params.get("extra", "")
    sjdb_overhang = snakemake.params.get("sjdbOverhang", "100")

    gtf = snakemake.input.get("gtf")
    if gtf is not None:
        gtf = "--sjdbGTFfile " + gtf
        sjdb_overhang = "--sjdbOverhang " + sjdb_overhang
    else:
        gtf = sjdb_overhang = ""

    makedirs(snakemake.output)

    shell(
        "STAR "  # Tool
        "--runMode genomeGenerate "  # Indexation mode
        "{extra} "  # Optional parameters
        "--runThreadN {snakemake.threads} "  # Number of threads
        "--genomeDir {snakemake.output} "  # Path to output
        "--genomeFastaFiles {snakemake.input.fasta} "  # Path to fasta files
        "{sjdb_overhang} "  # Read-len - 1
        "{gtf} "  # Highly recommended GTF
        "{log}"  # Logging
    )




* bio/star/align

.. code-block:: python

    __author__ = "Johannes Köster"
    __copyright__ = "Copyright 2016, Johannes Köster"
    __email__ = "koester@jimmy.harvard.edu"
    __license__ = "MIT"


    import os
    from snakemake.shell import shell

    extra = snakemake.params.get("extra", "")
    log = snakemake.log_fmt_shell(stdout=True, stderr=True)

    fq1 = snakemake.input.get("fq1")
    assert fq1 is not None, "input-> fq1 is a required input parameter"
    fq1 = (
        [snakemake.input.fq1]
        if isinstance(snakemake.input.fq1, str)
        else snakemake.input.fq1
    )
    fq2 = snakemake.input.get("fq2")
    if fq2:
        fq2 = (
            [snakemake.input.fq2]
            if isinstance(snakemake.input.fq2, str)
            else snakemake.input.fq2
        )
        assert len(fq1) == len(
            fq2
        ), "input-> equal number of files required for fq1 and fq2"
    input_str_fq1 = ",".join(fq1)
    input_str_fq2 = ",".join(fq2) if fq2 is not None else ""
    input_str = " ".join([input_str_fq1, input_str_fq2])

    if fq1[0].endswith(".gz"):
        readcmd = "--readFilesCommand zcat"
    else:
        readcmd = ""

    outprefix = os.path.dirname(snakemake.output[0]) + "/"

    shell(
        "STAR "
        "{extra} "
        "--runThreadN {snakemake.threads} "
        "--genomeDir {snakemake.params.index} "
        "--readFilesIn {input_str} "
        "{readcmd} "
        "--outFileNamePrefix {outprefix} "
        "--outStd Log "
        "{log}"
    )




* bio/arriba

.. code-block:: python

    __author__ = "Jan Forster"
    __copyright__ = "Copyright 2019, Jan Forster"
    __email__ = "j.forster@dkfz.de"
    __license__ = "MIT"


    import os
    from snakemake.shell import shell

    extra = snakemake.params.get("extra", "")
    log = snakemake.log_fmt_shell(stdout=True, stderr=True)

    discarded_fusions = snakemake.output.get("discarded", "")
    if discarded_fusions:
        discarded_cmd = "-O " + discarded_fusions
    else:
        discarded_cmd = ""

    blacklist = snakemake.params.get("blacklist")
    if blacklist:
        blacklist_cmd = "-b " + blacklist
    else:
        blacklist_cmd = ""

    known_fusions = snakemake.params.get("known_fusions")
    if known_fusions:
        known_cmd = "-k" + known_fusions
    else:
        known_cmd = ""

    sv_file = snakemake.params.get("sv_file")
    if sv_file:
        sv_cmd = "-d" + sv_file
    else:
        sv_cmd = ""

    shell(
        "arriba "
        "-x {snakemake.input.bam} "
        "-a {snakemake.input.genome} "
        "-g {snakemake.input.annotation} "
        "{blacklist_cmd} "
        "{known_cmd} "
        "{sv_cmd} "
        "-o {snakemake.output.fusions} "
        "{discarded_cmd} "
        "{extra} "
        "{log}"
    )




