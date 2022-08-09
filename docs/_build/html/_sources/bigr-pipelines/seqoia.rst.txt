.. _`seqOIA single-sample`:

SEQOIA SINGLE-SAMPLE
====================

Run seqOIA over single sample WES (tumor/normal) + RNA (tumor).

Usage
-----

In order to run the pipeline, use the following commands

.. code-block:: bash 

  # Go to your working directory

  cd /path/to/my/working/directory

  # Build a design file (see below)

  # Copy/paste the following line

  bash /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/seqoia/run.sh


Input/Output
------------


**Input:**

 
  
* Path to Fastq WES Tumor (R1 + R2)
  
 
  
* Path to Fastq WES Normal (R1 + R2)
  
 
  
* Path to Fastq RNA Tumor (R1 + R2)
  
 


**Output:**

 
  
* seqOIA results
  
 






Used wrappers
-------------

The following individual wrappers are used in this pipeline:


* :ref:`bio/BiGR/copy`


Please refer to each wrapper in above list for additional configuration parameters and information about the executed code.




Notes
-----

Prerequisites:

* A TSV formatted design file, *named 'design.tsv'* with the following columns:

.. list-table:: Desgin file format
  :widths: 33 33 33
  :header-rows: 1

  * - SAMPLE_ID
    - WES_T_UP
    - WES_T_DOWN
    - WES_N_UP
    - WES_N_DOWN
    - RNA_T_UP
    - RNA_T_DOWN
  * - Name of the Sample1
    - Path to upstream fastq WES Tumor
    - Path to downstream fastq WES Tumor
    - Path to upstream fastq WES Normal
    - Path to downstream fastq WES Normal
    - Path to upstream fastq RNA Tumor
    - Path to downstream fastq RNA Tumor
  * - Name of the Sample2
    - Path to upstream fastq WES Tumor
    - Path to downstream fastq WES Tumor
    - Path to upstream fastq WES Normal
    - Path to downstream fastq WES Normal
    - Path to upstream fastq RNA Tumor
    - Path to downstream fastq RNA Tumor
  * - ...
    - ...
    - ...
    - ...
    - ...
    - ...





Snakefile
---------

The pipeline contains the following steps:

.. code-block:: python

    #!/usr/bin/env python3
    # -*- coding: utf-8 -*-

    """This snakefile prepares each single sample run of SeqOIA and launch them"""

    from snakemake.utils import min_version
    from pathlib import Path
    from yaml import dump
    min_version("6.0")

    import sys

    worflow_source_dir = Path(snakemake.workflow.srcdir("."))
    common = str(worflow_source_dir / ".." / "common" / "python")
    sys.path.append(common)

    from file_manager import *
    from files_linker import *
    from write_yaml import *
    from messages import message

    logging.basicConfig(
        filename="snakemake.salmon_quant.log",
        filemode="w",
        level=logging.DEBUG
    )


    def link_seqoia_fq(
            sample_names: list[str],
            r1_paths: list[str],
            r2_paths: list[str],
            suffix: str
        ) -> dict[str, str]:
        """
        Case r2 are provided:
        Build a dictionnary containing the following pairs:
        original_r1_name: reads/{sample}.1.fq.gz
        original_r2_name: reads/{sample}.2.fq.gz

        Otherwise:
        Build a dictionnary containing the following fastq:
        original_name: reads/{sample}.fq.gz
        """
        link_dict = {}
        for sample, r1, r2 in zip(sample_names, r1_paths, r2_paths):
            link_dict[f"data_{sample}/A00000/22NNNN_ANNNNN_NNNN_AXXXXXXXXX/fastq/{sample}_{suffix}_S1_L1_11_R1_001.fastq.gz"] = r1
            link_dict[f"data_{sample}/A00000/22NNNN_ANNNNN_NNNN_AXXXXXXXXX/fastq/{sample}_{suffix}_S1_L1_11_R2_001.fastq.gz"] = r2
        return link_dict


    design = get_design(os.getcwd(), search_fastq_pairs)
    wes_t = link_seqoia_fq(
        design.SAMPLE_ID,
        design.WES_T_U,
        design.WES_T_D,
        suffix="WES-T"
    )

    wes_n = link_seqoia_fq(
        design.SAMPLE_ID,
        design.WES_N_U,
        design.WES_N_D,
        suffix="WGS-C"
    )

    rna_t = link_seqoia_fq(
        design.SAMPLE_ID,
        design.RNA_T_U,
        design.RNA_T_D,
        suffix="WTS"
    )


    sample_list = design.SAMPLE_ID.to_list()
    streams = ["1", "2"]
    datatypes = ["WTS", "WES-T", "WGS-C"]
    chromosome_list = list(map(str, range(1, 23))) + ["X", "Y", "MT"]

    wildcard_constraints:
        sample=r"|".join(sample_list),
        stream=r"|".join(streams),
        datatype=r"|".join(datatypes)


    rule target:
        input:
            expand(
                "{sample}.flag",
                #"data_{sample}/{sample}/{sample}_{datatype}_S1_L1_11",
                sample=sample_list,
                #datatype=datatypes
            )


    ##########################
    ### Build config files ###
    ##########################


    rule prepare_cluster_config:
        input:
            "/mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/seqoia/v.2.1.0.GRCh38.cluster.config.json"
        output:
            config = temp("v.2.1.0.GRCh38.cluster.config.json")
        message: "Gathering cluster configfile"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 128,
            time_min=lambda wildcards, attempt: attempt * 2,
            tmpdir="tmp"
        group:
            "General_config_files"
        log:
            "logs/prepare_cluster_config.log"
        params:
            "--verbose --checksum --human-readable --progress --partial"
        shell:
            "rsync {params} {input} {output} > {log} 2>&1"


    rule config_gleaves:
        output:
            ".config.gleaves"
        message:
            "Building config.gleaves file"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 128,
            time_min=lambda wildcards, attempt: attempt * 2,
            tmpdir="tmp"
        group:
            "General_config_files"
        shell:
            "touch {output}"


    rule prepare_sample_config:
        output:
            config = "v.2.1.0.GRCh38.{sample}.pipeline.config.json"
        message:
            "Building config file for {wildcards.sample}"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 512,
            time_min=lambda wildcards, attempt: attempt * 5,
            tmpdir="tmp"
        group: "sample"
        # conda:
        #     "env/seqoia.yaml"
        group:
            "sample"
        log:
            "logs/prepare_config/{sample}.log"
        script:
            "scripts/sample.config.py"


    #################################################
    ### Gather files from iRODS or mounting point ###
    #################################################


    rule gather_wes_tumor:
        output:
            temp("data_{sample}/A00000/22NNNN_ANNNNN_NNNN_AXXXXXXXXX/fastq/{sample}_WES-T_S1_L1_11_R{stream}_001.fastq.gz")
        message:
            "Gathering {wildcards.sample} fastq file (WES-T tumor {wildcards.stream})"
        threads: 1
        resources:
            mem_mb=lambda wildcard, attempt: min(attempt * 1024, 2048),
            time_min=lambda wildcard, attempt: attempt * 45,
            tmpdir="tmp"
        group: "sample"
        params:
            input=lambda wildcards, output: wes_t[output[0]],
            datatype="WES-T"
        log:
            "logs/bigr_copy/{sample}.WES-T.{stream}.log"
        wrapper:
            "bio/BiGR/copy"


    rule gather_wes_normal:
        output:
            temp("data_{sample}/A00000/22NNNN_ANNNNN_NNNN_AXXXXXXXXX/fastq/{sample}_WGS-C_S1_L1_11_R{stream}_001.fastq.gz")
        message:
            "Gathering {wildcards.sample} fastq file (WGS-C {wildcards.stream})"
        threads: 1
        resources:
            mem_mb=lambda wildcard, attempt: min(attempt * 1024, 2048),
            time_min=lambda wildcard, attempt: attempt * 45,
            tmpdir="tmp"
        group: "sample"
        params:
            input=lambda wildcards, output: wes_n[output[0]],
            datatype="WGS-C"
        log:
            "logs/bigr_copy/{sample}.WGS-C.{stream}.log"
        wrapper:
            "bio/BiGR/copy"


    rule gather_rna_tumoral:
        output:
            temp("data_{sample}/A00000/22NNNN_ANNNNN_NNNN_AXXXXXXXXX/fastq/{sample}_WTS_S1_L1_11_R{stream}_001.fastq.gz")
        message:
            "Gathering {wildcards.sample} fastq file (WTS tumor {wildcards.stream})"
        threads: 1
        resources:
            mem_mb=lambda wildcard, attempt: min(attempt * 1024, 2048),
            time_min=lambda wildcard, attempt: attempt * 45,
            tmpdir="tmp"
        group: "sample"
        params:
            input=lambda wildcards, output: rna_t[output[0]],
            datatype="WTS"
        log:
            "logs/bigr_copy/{sample}.WTS.{stream}.log"
        wrapper:
            "bio/BiGR/copy"


    ###########################
    ### Run SeqOIA pipeline ###
    ###########################


    rule create_empty_files:
        output:
            ped_genomiser=temp("data_{sample}/A00000/22NNNN_ANNNNN_NNNN_AXXXXXXXXX/ped_genomiser_generation.done"),
            sample_sheet=temp("data_{sample}/sample_sheet/22NNNN_ANNNNN_NNNN_AXXXXXXXXX.csv"),
            json=temp("data_{sample}/A00000/22NNNN_ANNNNN_NNNN_AXXXXXXXXX/json/{sample}_WGS-C.json"),
        message:
            "Building required empty files for {wildcards.sample}"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 128,
            time_min=lambda wildcards, attempt: attempt * 2,
            tmpdir="tmp"
        group: "sample"
        log:
            "logs/create_empty_files/{sample}.log"
        shell:
            "touch {output.ped_genomiser} {output.sample_sheet} {output.json} > {log} 2>&1"


    rule seqoia:
        input:
            sample_config = "v.2.1.0.GRCh38.{sample}.pipeline.config.json",
            pipeline_config = "v.2.1.0.GRCh38.cluster.config.json",
            fastq_paths = expand(
                "data_{sample}/A00000/22NNNN_ANNNNN_NNNN_AXXXXXXXXX/fastq/{sample}_{datatype}_S1_L1_11_R{stream}_001.fastq.gz",
                stream=streams,
                datatype=datatypes,
                sample="{sample}"
            ),
            snakefile="/mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/seqoia/scripts/Snakefile_analysis_nowgs_2.1",
            config_gleaves=".config.gleaves",
            ped_genomiser="data_{sample}/A00000/22NNNN_ANNNNN_NNNN_AXXXXXXXXX/ped_genomiser_generation.done",
            sample_sheet="data_{sample}/sample_sheet/22NNNN_ANNNNN_NNNN_AXXXXXXXXX.csv",
            json="data_{sample}/A00000/22NNNN_ANNNNN_NNNN_AXXXXXXXXX/json/{sample}_WGS-C.json",
        output:
            touch("{sample}.flag")
            #fusion_inspector=directory("data_{sample}/{sample}/{sample}_WTS_S1/star_fusion/fusion_inspector/annotator"),
            #logs=directory(expand(
            #    "data_{sample}/{sample}/log/{sample}_{datatype}_S1_L1_11",
            #    datatype=datatypes,
            #    allow_missing=True
            #)),
            #datatype_dir=directory(expand(
            #    "data_{sample}/{sample}/{sample}_{datatype}_S1_L1_11",
            #    datatype=datatypes,
            #    allow_missing=True
            #)),
            #charge_mut=directory("data_{sample}/{sample}/Charge_mutationnelle"),
            #chrom_dirs=directory(expand(
            #    "data_{sample}/{sample}/chr_{chromosome}",
            #    chromosome=chromosome_list,
            #    allow_missing=True
            #)),
            #fusion_dir=directory("data_{sample}/{sample}/Fusions"),
            #export_data_status="data_{sample}/{sample}/export_data_tar.done",
            #config_dir=directory("data_{sample}/{sample}/config"),
            #mut_sig=directory("data_{sample}/{sample}/Signature_mutationelle"),
            #sample_files=multiext(
            #    "data_{sample}/{sample}/{sample}",
            #    "_final_temp.vcf",
            #    "_final.vcf",
            #    "_importation_gleaves.done",
            #    ".json",
            #    "_name.vcf",
            #    "_qc.csv",
            #    "_qc_lab.csv",
            #    "_qc_lab_temp.csv",
            #    "_qc_validation.done",
            #    ".tar.gz"
            #)
        message:
            "Running SeqOIA single sample no WGS on {wildcards.sample}"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 1024 * 2,
            time_min=lambda wildcards, attempt: attempt * 60 * 10,
            # mem_mb = 4096,
            # time_min = 5,
            tmpdir="tmp"
        shadow: "minimal"
        params:
            mk = "--parents --verbose",
            ln = "--symbolic --force --relative",
            snake_args = "--jobs 30 --rerun-incomplete --keep-going --max-jobs-per-second 1",
            sbatch = "'sbatch -p {cluster.partition} -t {cluster.time} --cpus-per-task {cluster.cpu} --mem {cluster.mem} --exclude=n01,n02,n03,n04,n05,n06'",
            fusion_inspector=directory("data_{sample}/{sample}/{sample}_WTS_S1/star_fusion/fusion_inspector/annotator"),
            logs=directory(expand(
                "data_{sample}/{sample}/log/{sample}_{datatype}_S1_L1_11",
                datatype=datatypes,
                allow_missing=True
            )),
            datatype_dir=directory(expand(
                "data_{sample}/{sample}/{sample}_{datatype}_S1_L1_11",
                datatype=datatypes,
                allow_missing=True
            )),
        conda:
            "env/seqoia.yaml"
        log:
            "logs/seqOIA/{sample}.log"
        shell:
            "mkdir {params.mk} {params.fusion_inspector} {params.logs} .snakemake "
            #"mkdir {params.mk} {output.fusion_inspector} {output.logs} .snakemake "
            "{params.datatype_dir} data_{wildcards.sample}/sample_sheet/"
            #"{output.datatype_dir} data_{wildcards.sample}/sample_sheet/"
            "data_{wildcards.sample}/json/ > {log} 2>&1 && "
            "tree -a >> {log} 2>&1 && "
            "pwd >> {log} 2>&1 && "
            "snakemake --snakefile {input.snakefile} "
            "--configfile {input.sample_config} "
            "--cluster-config {input.pipeline_config} "
            "--cluster {params.sbatch} "
            "{params.snake_args} --unlock >> {log} 2>&1 && "
            "snakemake --snakefile {input.snakefile} "
            "--configfile {input.sample_config} "
            "--cluster-config {input.pipeline_config} "
            "--cluster {params.sbatch} "
            "{params.snake_args} >> {log} 2>&1 "




Authors
-------


* Thibault Dayris

* Marc Deloger
