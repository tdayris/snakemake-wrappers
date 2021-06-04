.. _`EaCoN Post Process`:

EACON POST PROCESS
==================

All EaCoN operations following the Oncoscan or Cytoscan Process


Example
-------

This meta-wrapper can be used by integrating the following into your workflow:

.. code-block:: python

    # This snakefile should be called through other pipelines !
    # Use either :
    # * cytoscan_eacon
    # * oncoscan_eacon
    # WARNING: Config is supposed to be filled !

    def get_gamma_files(config):
        eval = [
            "{sample}/{segmenter}/ASCN/{sample}.gammaEval.{ext}".format(
                sample="{sample}", segmenter=config["params"]["segmenter"], ext=ext
            )
            for ext in ["png", "txt"]
        ]

        # Sometimes EaCoN formats floats with trailing zeros, sometimes it does not
        gamma_formats = zip(
            [f"{x/100:.2f}" for x in range(35, 95, 5)],
            [f"{x/100}" for x in range(35, 95, 5)]
        )

        results = [
            "{sample}/{segmenter}/ASCN/gamma{g1}/{sample}.gamma{g2}{files}".format(
                sample="{sample}",
                segmenter=config["params"]["segmenter"],
                g1=g1,
                g2=g2,
                files=ext
            )
            for g1, g2 in gamma_formats
            for ext in [
                ".cn",
                "_model.txt",
            ]
        ]
        results += [
            "{sample}/{segmenter}/ASCN/gamma{g1}/{sample}.{files}".format(
                sample="{sample}",
                segmenter=config["params"]["segmenter"],
                g1=g1,
                g2=g2,
                files=ext
            )
            for g1, g2 in gamma_formats
            for ext in [
                ".rawprofile.png",
                ".Rorschach.clown.png",
                ".ASCN.ASCAT.png",
                ".ASCN.ASCAT.RDS",
                ".TCNvsL2R.png"
            ]
        ]

        return {
            "eval": eval,
            "gamma_results": results
        }

    rule eacon_genomic_instability:
        input:
            **get_gamma_files(config=config)
        output:
            "{sample}/{sample}_GIS_from_best_gamma.txt"
        threads: 1
        resources:
            time_min=lambda wildcards, attempt: attempt * 35,
            mem_mb=lambda wildcards, attempt: attempt * 2 * 1024
        log:
            "logs/EaCoN/{sample}/genomic_instability.log"
        params:
            indir = "{sample}",  # Required!
            genome = config["params"]["genome"],
            segmenter = config["params"]["segmenter"]
        wrapper:
            "/bio/eacon/instability"

    rule eacon_annotate:
        input:
            rds = "{sample}/ASCAT/L2R/{sample}.SEG.ASCAT.RDS",
            grd = "grd",
            ldb = "databases"
        output:
            multiext(
                "{sample}/ASCAT/L2R/{sample}",
                ".BAF.png",
                ".Cut.acbs",
                ".Instab.txt",
                ".INT.png",
                ".L2R.G.png",
                ".L2R.K.png",
                ".TargetGenes.txt",
                ".TruncatedGenes.txt"
            ),
            expand(
                "{sample}/ASCAT/L2R/chromosomes/{chromosome}",
                chromosome = [
                    f"chr{i}.png" for i in list(map(str, range(1, 23))) + ["X", "Y"]
                ],
                allow_missing=True
            ),
            "{sample}/ASCAT/L2R/{sample}.REPORT.html",
            directory("{sample}/ASCAT/L2R/{sample}_solo.hg19")
        threads: 1
        resources:
            time_min=lambda wildcards, attempt: attempt * 35,
            mem_mb=lambda wildcards, attempt: attempt * 3 * 1024
        log:
            "logs/EaCoN/{sample}/annotate.log"
        wrapper:
            "/bio/eacon/annotate"


    rule eacon_databases:
        output:
            databases = directory("databases")
        log:
            "logs/EaCoN/databases.log"
        cache: True
        wrapper:
            "/bio/eacon/databases"


    rule eacon_grd:
        output:
            grd = "grd"
        log:
            "logs/EaCoN/grd.log"
        cache: True
        wrapper:
            "/bio/eacon/databases"


    rule eacon_ascn:
        input:
            rds = "{sample}/ASCAT/L2R/{sample}.SEG.ASCAT.RDS"
        output:
            **get_gamma_files(config=config)
        threads: 1
        resources:
            time_min=lambda wildcards, attempt: attempt * 35,
            mem_mb=lambda wildcards, attempt: attempt * 3 * 1024
        log:
            "logs/EaCoN/{sample}/ascn.log"
        params:
            extra = config["extra"]["EaCoN_ascn"]
        wrapper:
            "/bio/eacon/ascn"

    rule EaCoN_segment:
        input:
            rds = "{sample}/{sample}_{arraytype}_{genome}_processed.RDS".format(
                sample="{sample}",
                nar=config["params"]["nar"],
                genome=config["params"]["genome"],
                arraytype = config["params"]["arraytype"]
            )
        output:
            files = multiext(
                "{sample}/ASCAT/L2R/{sample}",
                ".Cut.cbs",
                ".NoCut.cbs",
                ".Rorschach.png",
                ".SegmentedBAF.txt"
            ),
            rds = "{sample}/ASCAT/L2R/{sample}.SEG.ASCAT.RDS",
            png = "{sample}/ASCAT/L2R/{sample}.SEG.ASCAT.png"
        threads: 1
        resources:
            time_min=lambda wildcards, attempt: attempt * 35,
            mem_mb=lambda wildcards, attempt: attempt * 3 * 1024
        log:
            "logs/EaCoN/{sample}/segment.log"
        params:
            extra = config["extra"]["EaCoN_segment"]
        wrapper:
            "/bio/eacon/segment"

Note that input, output and log file paths can be chosen freely, as long as the dependencies between the rules remain as listed here.
For additional parameters in each individual wrapper, please refer to their corresponding documentation (see links below).

When running with

.. code-block:: bash

    snakemake --use-conda

the software dependencies will be automatically deployed into an isolated environment before execution.



Used wrappers
---------------------

The following individual wrappers are used in this meta-wrapper:


* :ref:`master/bio/eacon/segment`

* :ref:`master/bio/eacon/databases`

* :ref:`master/bio/eacon/ascn`

* :ref:`master/bio/eacon/annotate`

* :ref:`master/bio/eacon/instability`


Please refer to each wrapper in above list for additional configuration parameters and information about the executed code.






Notes
-----

This meta wrapper is called through the EaCoN_OncoScan or the EaCoN_Cytoscan
meta-wrappers

Configuration architecture is duscissed in these meta-wrappers' pages.




Authors
-------


* Thibault Dayris

* Bastien Job

