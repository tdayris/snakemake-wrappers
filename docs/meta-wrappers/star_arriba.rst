.. _`star-arriba`:

STAR-ARRIBA
===========

A subworkflow for fusion detection from RNA-seq data with ``arriba``. The fusion calling is based on splice-aware, chimeric alignments done with ``STAR``. ``STAR`` is used with specific parameters to ensure optimal functionality of the ``arriba`` fusion detection, for details, see the `documentation <https://arriba.readthedocs.io/en/latest/workflow/>`_.



Example
-------

This meta-wrapper can be used by integrating the following into your workflow:

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
            "0.72.0-501-g28151774c/bio/star/index"

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
            "0.72.0-501-g28151774c/bio/star/align"

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
            "0.72.0-501-g28151774c/bio/arriba"

Note that input, output and log file paths can be chosen freely, as long as the dependencies between the rules remain as listed here.
For additional parameters in each individual wrapper, please refer to their corresponding documentation (see links below).

When running with

.. code-block:: bash

    snakemake --use-conda

the software dependencies will be automatically deployed into an isolated environment before execution.



Used wrappers
---------------------

The following individual wrappers are used in this meta-wrapper:


* :ref:`bio/star/index`

* :ref:`bio/star/align`

* :ref:`bio/arriba`


Please refer to each wrapper in above list for additional configuration parameters and information about the executed code.







Authors
-------


* Jan Forster

