#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""This snakefile prepares each single sample run of SeqOIA and launch them"""

from snakemake.utils import min_version
from pathlib import Path
from yaml import dump
from json import load
min_version("6.0")

import sys

worflow_source_dir = Path(next(iter(workflow.get_sources()))).absolute().parent
common = str(worflow_source_dir / "../common/python")
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
            "data_{sample}/{sample}/{sample}_{datatype}_S1_L1_11",
            sample=sample_list,
            datatype=datatypes
        )


##########################
### Build config files ###
##########################

#
# rule prepare_cluster_config:
#     input:
#         "/mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/seqoia/v.2.1.0.GRCh38.cluster.config.json"
#     output:
#         config = temp("v.2.1.0.GRCh38.cluster.config.json")
#     message: "Gathering cluster configfile"
#     threads: 1
#     resources:
#         mem_mb=lambda wildcards, attempt: attempt * 128,
#         time_min=lambda wildcards, attempt: attempt * 2,
#         tmpdir="tmp"
#     group:
#         "General_config_files"
#     log:
#         "logs/prepare_cluster_config.log"
#     params:
#         "--verbose --checksum --human-readable --progress --partial"
#     shell:
#         "rsync {params} {input} {output} > {log} 2>&1"


# rule config_gleaves:
#     output:
#         ".config.gleaves"
#     message:
#         "Building config.gleaves file"
#     threads: 1
#     resources:
#         mem_mb=lambda wildcards, attempt: attempt * 128,
#         time_min=lambda wildcards, attempt: attempt * 2,
#         tmpdir="tmp"
#     group:
#         "General_config_files"
#     shell:
#         "touch {output}"


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
        "Sample_Preparation"
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
        "Gathering {wildcards.sample} fastq file (WES tumor {wildcards.stream})"
    threads: 1
    resources:
        mem_mb=lambda wildcard, attempt: min(attempt * 1024, 2048),
        time_min=lambda wildcard, attempt: attempt * 45
    group: "sample"
    params:
        input=lambda wildcards, output: wes_t[output[0]]
    log:
        "logs/bigr_copy/{sample}.{stream}.log"
    wrapper:
        "bio/BiGR/copy"


rule gather_wes_normal:
    output:
        temp("data_{sample}/A00000/22NNNN_ANNNNN_NNNN_AXXXXXXXXX/fastq/{sample}_WGS-C_S1_L1_11_R{stream}_001.fastq.gz")
    message:
        "Gathering {wildcards.sample} fastq file (WES normal {wildcards.stream})"
    threads: 1
    resources:
        mem_mb=lambda wildcard, attempt: min(attempt * 1024, 2048),
        time_min=lambda wildcard, attempt: attempt * 45
    group: "sample"
    params:
        input=lambda wildcards, output: wes_n[output[0]]
    log:
        "logs/bigr_copy/{sample}.{stream}.log"
    wrapper:
        "bio/BiGR/copy"


rule gather_rna_tumoral:
    output:
        temp("data_{sample}/A00000/22NNNN_ANNNNN_NNNN_AXXXXXXXXX/fastq/{sample}_WTS_S1_L1_11_R{stream}_001.fastq.gz")
    message:
        "Gathering {wildcards.sample} fastq file (RNA tumor {wildcards.stream})"
    threads: 1
    resources:
        mem_mb=lambda wildcard, attempt: min(attempt * 1024, 2048),
        time_min=lambda wildcard, attempt: attempt * 45
    group: "sample"
    params:
        input=lambda wildcards, output: rna_t[output[0]]
    log:
        "logs/bigr_copy/{sample}.{stream}.log"
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


with open("v.2.1.0.GRCh38.cluster.config.json") as seqoyaml:
    seqoia_config = load(seqoyaml)


module seqoia_workflow:
    snakefile: "/mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/seqoia/scripts/Snakefile_analysis_nowgs_2.1"
    config: seqoia_config


use rule * from seqoia_workflow
#
# rule seqoia:
#     input:
#         sample_config = "v.2.1.0.GRCh38.{sample}.pipeline.config.json",
#         pipeline_config = "v.2.1.0.GRCh38.cluster.config.json",
#         fastq_paths = expand(
#             "data_{sample}/A00000/22NNNN_ANNNNN_NNNN_AXXXXXXXXX/fastq/{sample}_{datatype}_S1_L1_11_R{stream}_001.fastq.gz",
#             stream=streams,
#             datatype=datatypes,
#             sample="{sample}"
#         ),
#         snakefile="/mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/seqoia/scripts/Snakefile_analysis_nowgs_2.1",
#         config_gleaves=".config.gleaves",
#         ped_genomiser="data_{sample}/A00000/22NNNN_ANNNNN_NNNN_AXXXXXXXXX/ped_genomiser_generation.done",
#         sample_sheet="data_{sample}/sample_sheet/22NNNN_ANNNNN_NNNN_AXXXXXXXXX.csv",
#         json="data_{sample}/A00000/22NNNN_ANNNNN_NNNN_AXXXXXXXXX/json/{sample}_WGS-C.json",
#     output:
#         fusion_inspector=directory("data_{sample}/{sample}/{sample}_WTS_S1/star_fusion/fusion_inspector/annotator"),
#         logs=directory(expand(
#             "data_{sample}/{sample}/log/{sample}_{datatype}_S1_L1_11",
#             datatype=datatypes,
#             allow_missing=True
#         )),
#         datatype_dir=directory(expand(
#             "data_{sample}/{sample}/{sample}_{datatype}_S1_L1_11",
#             datatype=datatypes,
#             allow_missing=True
#         )),
#         charge_mut=directory("data_{sample}/{sample}/Charge_mutationnelle"),
#         chrom_dirs=directory(expand(
#             "data_{sample}/{sample}/chr_{chromosome}",
#             chromosome=chromosome_list,
#             allow_missing=True
#         )),
#         fusion_dir=directory("data_{sample}/{sample}/Fusions"),
#         export_data_status="data_{sample}/{sample}/export_data_tar.done",
#         config_dir=directory("data_{sample}/{sample}/config"),
#         mut_sig=directory("data_{sample}/{sample}/Signature_mutationelle"),
#         sample_files=multiext(
#             "data_{sample}/{sample}/{sample}",
#             "_final_temp.vcf",
#             "_final.vcf",
#             "_importation_gleaves.done",
#             ".json",
#             "_name.vcf",
#             "_qc.csv",
#             "_qc_lab.csv",
#             "_qc_lab_temp.csv",
#             "_qc_validation.done",
#             ".tar.gz"
#         )
#     message:
#         "Running SeqOIA single sample no WGS on {wildcards.sample}"
#     threads: 1
#     resources:
#         mem_mb=lambda wildcards, attempt: attempt * 1024 * 2,
#         time_min=lambda wildcards, attempt: attempt * 60 * 10,
#         # mem_mb = 4096,
#         # time_min = 5,
#         tmpdir="tmp"
#     shadow: "minimal"
#     group: "sample"
#     params:
#         mk = "--parents --verbose",
#         ln = "--symbolic --force --relative",
#         snake_args = "--jobs 30 --rerun-incomplete --keep-going --max-jobs-per-second 1",
#         sbatch = "'sbatch -p {cluster.partition} -t {cluster.time} --cpus-per-task {cluster.cpu} --mem {cluster.mem} --exclude=n01,n02,n03,n04,n05,n06'"
#     conda:
#         "env/seqoia.yaml"
#     log:
#         "logs/seqOIA/{sample}.log"
#     shell:
#         "mkdir {params.mk} {output.fusion_inspector} {output.logs} .snakemake"
#         "{output.datatype_dir} data_{wildcards.sample}/sample_sheet/"
#         "data_{wildcards.sample}/json/ > {log} 2>&1 && "
#         "tree -a >> {log} 2>&1 && "
#         "pwd >> {log} 2>&1 && "
#         "snakemake --snakefile {input.snakefile} "
#         "--configfile {input.sample_config} "
#         "--cluster-config {input.pipeline_config} "
#         "--cluster {params.sbatch} "
#         "{params.snake_args} --unlock >> {log} 2>&1 && "
#         "snakemake --snakefile {input.snakefile} "
#         "--configfile {input.sample_config} "
#         "--cluster-config {input.pipeline_config} "
#         "--cluster {params.sbatch} "
#         "{params.snake_args} >> {log} 2>&1 "
