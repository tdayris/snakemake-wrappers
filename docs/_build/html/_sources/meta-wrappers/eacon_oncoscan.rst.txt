.. _`EaCoN Cytoscan`:

EACON CYTOSCAN
==============

Analyse CNV on CytoScanHD Arrays with EaCoN


Example
-------

This meta-wrapper can be used by integrating the following into your workflow:

.. code-block:: python

    from snakemake.utils import min_version
    from pathlib import Path
    from yaml import dump
    min_version("6.0")


    # If not config file is provided, a default configfile will be created.
    default_oncoscan_process_config = {
        "samples": [""],
        "params": {
            "arraytype": "OncoScan_CNV",
            "genome": "hg19",
            "nar": "na33.r2",
            "segmenter": "ASCAT"
        },
        "extra": {
            "EaCoN_process": 'apt.build = "na33.r2", force = TRUE',
            "EaCoN_segment": (
                'segmenter = "ASCAT", smooth.k = NULL, BAF.filter = 0.9, '
                'SER.pen = 40, nrf = 0.5, penalty = 50, force = TRUE'
            ),
            "EaCoN_ascn": 'force = TRUE',
            "EaCoN_gis": ""
        }
    }

    config_path = Path("config.yaml")
    if not config_path.exists():
        # Search for CEL files in order to get sample names
        sample_identifiers = [
            cel.name[:-len("_A.CEL")]
            for cel in config_path.parent.iterdir()
            if cel.name.endswith("_A.CEL")
        ]
        # Update config file with sample names
        default_oncoscan_process_config["samples"] = sample_identifiers

        # Save config file
        with config_path.open("w") as config_stream:
            config_stream.write(
                dump(default_oncoscan_process_config, default_flow_style=False)
            )

    # The configfile name is hard coded. It should be: "config.yaml"
    configfile: str(config_path)

    module post_process_eacon:
        snakefile: "../../eacon_post_process/test/Snakefile"
        config: config

    print(config)

    rule default_oncoscan_process_all:
        input:
            # EaCoN models
            ascn = expand(
                os.sep.join(["{sample}", config["params"]["segmenter"],
                             "ASCN", "{sample}.gammaEval.png"]),
                sample=config["samples"]
            ),
            # EaCoN annotate
            html = expand(
                os.path.sep.join([
                    "{sample}", config["params"]["segmenter"], "L2R",
                    "{sample}.REPORT.html"
                ]),
                sample=config["samples"]
            ),
            # EaCoN new instability scoring feature
            #instability = expand(
            #    "{sample}/{sample}_GIS_from_best_gamma.txt",
            #    sample=config["samples"]
            #)

    # Import all rules from the eacon_post_process meta wrapper
    use rule * from post_process_eacon as post_process_*

    rule eacon_oncoscan_process:
        input:
            install = "sources",
            ATChannelCel = "{sample}_A.CEL",
            GCChannelCel = "{sample}_C.CEL"
        output:
            qc_txt = "{sample}/{sample}_2.4.0_{nar}.qc.txt".format(
                sample="{sample}", nar=config["params"]["nar"]
            ),
            log = "{sample}/{sample}_2.4.0_{nar}.log".format(
                sample="{sample}", nar=config["params"]["nar"]
            ),
            txt = "{sample}/{sample}_pairs.txt",
            png = "{sample}/{sample}_{arraytype}_{genome}_rawplot.png".format(
                sample="{sample}",
                arraytype=config["params"]["arraytype"],
                genome=config["params"]["genome"]
            ),
            rds = "{sample}/{sample}_{arraytype}_{genome}_processed.RDS".format(
                sample="{sample}",
                genome=config["params"]["genome"],
                arraytype = config["params"]["arraytype"]
            ),
            pairs = "{sample}/{sample}_2.4.0_{nar}.paircheck.txt".format(
                sample="{sample}", nar=config["params"]["nar"]
            )
        threads: 1
        resources:
            time_min=lambda wildcards, attempt: attempt * 50,
            mem_mb=lambda wildcards, attempt: attempt * 5 * 1024
        params:
            extra = config["extra"]["EaCoN_process"]
        log:
            "logs/EaCoN/{sample}/oncoscan_process.log"
        wrapper:
            "/bio/eacon/oncoscan_process"


    rule eacon_install:
        output:
            directory("sources")
        params:
            OncoScan = True,
            OncoScanCNV = True,
            CytoScan750K = True,
            CytoScanHD = True,
            genomewide = False,  # WARNING: Genome wide information not installed
            norm = True,
            EaCoN_dev = True,
            EaCoN_chromosomes = True
        cache: True
        threads: 1
        resources:
            time_min = lambda wildcards, attempt: attempt * 480,
            mem_mb = lambda wildcards, attempt: attempt * 4096
        log:
            "logs/EaCoN/install.log"
        wrapper:
            "/bio/eacon/install"

Note that input, output and log file paths can be chosen freely, as long as the dependencies between the rules remain as listed here.
For additional parameters in each individual wrapper, please refer to their corresponding documentation (see links below).

When running with

.. code-block:: bash

    snakemake --use-conda

the software dependencies will be automatically deployed into an isolated environment before execution.



Used wrappers
---------------------

The following individual wrappers are used in this meta-wrapper:


* :ref:`master/bio/eacon/install`

* :ref:`master/bio/eacon/cytoscan_process`

* :ref:`master/bio/eacon/segment`

* :ref:`master/bio/eacon/databases`

* :ref:`master/bio/eacon/ascn`

* :ref:`master/bio/eacon/annotate`

* :ref:`master/bio/eacon/instability`


Please refer to each wrapper in above list for additional configuration parameters and information about the executed code.






Notes
-----

notes




Authors
-------


* Thibault Dayris

* Bastien Job

