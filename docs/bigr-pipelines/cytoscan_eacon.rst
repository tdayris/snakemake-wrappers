.. _`EaCoN_Cytoscan`:

EACON_CYTOSCAN
==============

Analyse Cytoscans with EaCoN on Flamingo

Usage
-----

In order to run the pipeline, use the following commands

.. code-block:: bash 

  # Go to your directory containing CEL files

  cd /path/to/CEL/dir

  # Copy/paster the following line

  bash /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/cytoscan_eacon/run.sh


Input/Output
------------


**Input:**

 
  
* CEL files
  
 


**Output:**

 
  
* Complete EaCoN anaysis directory
  
 
  
* HTML report
  
 





Used meta-wrappers
------------------

The following individual meta-wrappers are used in this pipeline:


* :ref:`meta/bio/eacon_cytoscan`

* :ref:`meta/bio/eacon_post_process`


Please refer to each meta-wrapper in above list for additional configuration parameters and information about the executed code.







Snakefile
---------

The pipeline contains the following steps:

.. code-block:: python

    from snakemake.utils import min_version
    from pathlib import Path
    from yaml import dump
    min_version("6.0")


    # If not config file is provided, a default configfile will be created.
    default_cytoscan_process_config = {
        "samples": [""],
        "params": {
            "arraytype": "CytoScanHD_Array",
            "genome": "hg19",
            "nar": "na33.r4",
            "segmenter": "ASCAT"
        },
        "extra": {
            "EaCoN_process": 'apt.build = "na33.r4", force = TRUE',
            "EaCoN_segment": (
                'segmenter = "ASCAT", smooth.k = 5, BAF.filter = 0.75, '
                'SER.pen = 20, nrf = 1.0, penalty = 50, force = TRUE'
            ),
            "EaCoN_ascn": 'force = TRUE',
            "EaCoN_gis": ""
        }
    }

    config_path = Path("config.yaml")
    if not config_path.exists():
        # Search for CEL files in order to get sample names
        sample_identifiers = [
            cel.name[:-len(".CEL")]
            for cel in config_path.parent.iterdir()
            if cel.name.endswith(".CEL")
        ]
        # Update config file with sample names
        default_cytoscan_process_config["samples"] = sample_identifiers

        # Save config file
        with config_path.open("w") as config_stream:
            config_stream.write(
                dump(default_cytoscan_process_config, default_flow_style=False)
            )

    # The configfile name is hard coded. It should be: "config.yaml"
    configfile: str(config_path)

    module post_process_eacon:
        snakefile: "../../meta/bio/eacon_post_process/test/Snakefile"
        config: config

    ruleorder: eacon_annotate > post_process_eacon_annotate

    rule default_cytoscan_process_all:
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

    use rule eacon_annotate from post_process_eacon with:
        input:
            rds = "{sample}/ASCAT/L2R/{sample}.SEG.ASCAT.RDS",
            grd = "/mnt/beegfs/database/bioinfo/Index_DB/EaCoN/scripts/grd",
            ldb = "/mnt/beegfs/database/bioinfo/Index_DB/EaCoN/databases"

    rule eacon_cytoscan_process:
        input:
            install = 'install.ok',
            cel = "{sample}.CEL",
        output:
            qc_txt = "{sample}/{sample}_2.4.0_{nar}.qc.txt".format(
                sample="{sample}",
                nar=config["params"]["nar"]
            ),
            log = "{sample}/{sample}_2.4.0_{nar}.log".format(
                sample="{sample}",
                nar=config["params"]["nar"]
            ),
            txt = "{sample}/{sample}_CELfile.txt",
            png = "{sample}/{sample}_{arraytype}_{genome}_rawplot.png".format(
                sample="{sample}",
                arraytype=config["params"]["arraytype"],
                genome=config["params"]["genome"]
            ),
            rds = "{sample}/{sample}_{arraytype}_{genome}_processed.RDS".format(
                sample="{sample}",
                genome=config["params"]["genome"],
                arraytype = config["params"]["arraytype"]
            )
        threads: 1
        resources:
            time_min=lambda wildcards, attempt: attempt * 50,
            mem_mb=lambda wildcards, attempt: attempt * 5 * 1024
        params:
            extra = config["extra"]["EaCoN_process"]
        log:
            "logs/EaCoN/{sample}/cytoscan_process.log"
        wrapper:
            "/bio/eacon/cytoscan_process"


    rule eacon_install:
        input:
            r_packages = [
                "/mnt/beegfs/database/bioinfo/Index_DB/EaCoN/packages/affy.CN.norm.data_0.1.2.tar.gz",
                "/mnt/beegfs/database/bioinfo/Index_DB/EaCoN/packages/CytoScan750K.Array.na33.r4_0.1.0.tar.gz",
                "/mnt/beegfs/database/bioinfo/Index_DB/EaCoN/packages/CytoScan750K.Array.na36.r1_0.1.0.tar.gz",
                "/mnt/beegfs/database/bioinfo/Index_DB/EaCoN/packages/CytoScanHD.Array.na33.r4_0.1.0.tar.gz",
                "/mnt/beegfs/database/bioinfo/Index_DB/EaCoN/packages/CytoScanHD.Array.na36.r1_0.1.0.tar.gz",
                "/mnt/beegfs/database/bioinfo/Index_DB/EaCoN/packages/OncoScanCNV.na33.r2_0.1.0.tar.gz",
                "/mnt/beegfs/database/bioinfo/Index_DB/EaCoN/packages/OncoScanCNV.na36.r1_0.1.0.tar.gz",
                "/mnt/beegfs/database/bioinfo/Index_DB/EaCoN/packages/OncoScan.na33.r4_0.1.0.tar.gz",
                "/mnt/beegfs/database/bioinfo/Index_DB/EaCoN/packages/OncoScan.na36.r1_0.1.0.tar.gz",
                "/mnt/beegfs/database/bioinfo/Index_DB/EaCoN/packages/rcnorm_0.1.5.tar.gz"
            ],
            git_packages = [
                "/mnt/beegfs/database/bioinfo/Index_DB/EaCoN/packages/EaCoN_0.75.0-772-g3f7df6e90.tar.gz",
                "/mnt/beegfs/database/bioinfo/Index_DB/EaCoN/packages/EaCoN_Chromosomes.tar.gz",
                "/mnt/beegfs/database/bioinfo/Index_DB/EaCoN/packages/apt.cytoscan.2.4.0.tar.gz",
                "/mnt/beegfs/database/bioinfo/Index_DB/EaCoN/packages/apt.oncoscan.2.4.0.tar.gz"
            ]
        output:
            temp(touch('install.ok'))
        cache: True
        threads: 1
        resources:
            time_min = lambda wildcards, attempt: attempt * 480,
            mem_mb = lambda wildcards, attempt: attempt * 4096
        log:
            "logs/EaCoN/install.log"
        wrapper:
            "/bio/eacon/install_local"




Authors
-------


* Thibault Dayris

* Bastien Job

* Gérôme Jules-Clément
