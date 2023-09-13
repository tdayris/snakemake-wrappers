"""
This snakefile contains python functions and globals
Skip if you're looking for rules
"""

#####################################
# Check Snakemake version           #
# Import search-and-buil functions  #
# for config and design             #
#####################################

# Official libraries
import os
import functools
import itertools

from snakemake.utils import min_version
from pathlib import Path
from yaml import dump

min_version("7.5")

import sys
import pandas

# My own libraries
workflow_source_dir = Path(snakemake.workflow.srcdir(".."))
common = str(workflow_source_dir / ".." / "common" / "python")
sys.path.append(common)

from dataframes import *
from file_manager import *
from files_linker import *
from graphics import *
from write_yaml import *
from reservation import *
from messages import message

#####################
# Setup environment #
#####################

# Save output stream in a file
logging.basicConfig(filename="snakemake.CTC_Integragen_Legacy.log", filemode="w", level=logging.DEBUG)
logging.info("Additional utils loaded")


# Find and load configfile
default_config = read_yaml(workflow_source_dir / "config.hg38.yaml")


configfile: get_config(default_config=default_config)


logging.info("Config file loaded")


# Load design file and duplicate sample id as row name
design = get_design(dirpath=os.getcwd(), search_func=search_mapping)
design.set_index("Sample_id", inplace=True)
design["Sample_id"] = design.index.tolist()
design = design.astype("str")
logging.info("Design file loaded")

##################################
# Setup globals and fix wildcards #
##################################

# Links fastq paths provided by users and fastq paths used in this pipeline
# this is done in order to handle iRODS paths.
logging.info("Building globals...")

# def parse_design(
#     design: pandas.DataFrame, prefix: str = "data_input", suffix: str = "bam"
# ) -> dict[str, str]:
#     "We assume Baseline and WBC correspondig to a sample always come first."

#     logging.info("Parsing design...")
#     link_bams = {}
#     sample_list = []
#     baseline_sample_list = []
#     wbc_sample_list = []
#     link_sample_baseline = {}
#     last_baseline = None
#     last_wbc = None

#     row_iter = iter(design.iterrows())
#     row = next(row_iter, None)

#     while row is not None:
#         row = row[1]
#         sample = row["Sample_id"]
#         if row["Status"].lower() == "baseline":
#             link_bams[f"{prefix}/{sample}.baseline.{suffix}"] = row["bam"]
#             logging.debug(f"New baseline added for {sample} in general")
#             baseline_sample_list.append(sample)
#             last_baseline = f"{sample}.baseline"

#         elif row["Status"].lower() == "wbc":
#             manip = row["Manip"]
#             kit = row["Version"]
#             sample_id = f"{sample}_{kit}_M{manip}"
#             wbc_sample_list.append(sample_id)
#             last_wbc = f"{sample_id}.wbc"

#             link_bams[f"{prefix}/{sample_id}.wbc.{suffix}"] = row["bam"]
#             logging.debug(
#                 f"New WBC added for {sample_id} precisely, all replicates concerned."
#             )

#         elif row["Status"].lower() == "ctc":
#             manip = row["Manip"]
#             kit = row["Version"]
#             replicate = row["Replicate"]
#             raw_sample_id = f"{sample}_{kit}_M{manip}"
#             sample_id = f"{raw_sample_id}_{replicate}"

#             sample_list.append(sample_id)
#             link_bams[f"{prefix}/{sample_id}.ctc.{suffix}"] = row["bam"]
#             link_sample_baseline[sample_id] = {
#                 "ctc": f"sambamba/markdup/{sample_id}.ctc.{suffix}",
#                 "wbc": f"sambamba/markdup/{last_wbc}.{suffix}",
#                 "baseline": f"sambamba/markdup/{last_baseline}.{suffix}",
#             }
#             logging.debug(
#                 f"New CTC added {raw_sample_id}, replicate number {replicate}."
#             )

#         else:
#             raise ValueError(row)

#         row = next(row_iter, None)


#     return link_bams, sample_list, link_sample_baseline, baseline_sample_list, wbc_sample_list

# def get_baseline(wildcards):
#     return link_sample_baseline[wildcards.sample]["baseline"]

def get_baseline(wildcards):
    # Optimal case.
    if "nb" in wildcards.keys():
        if "manip" in wildcards.keys():
            if "version" in wildcards.keys():
                if "sample" in wildcards.keys():
                    return design[
                        (design["Sample_id"] == str(wildcards.sample)) & 
                        (design["Version"] == str(wildcards.version)) &
                        (design["Manip"] == str(wildcards.manip)) &
                        (design["NB"] == str(wildcards.nb))
                    ].Baseline.to_list()[0]

    # There is only one baseline per sample, so
    # missing "nb" should not be a problem.
    if "manip" in wildcards.keys():
        if "version" in wildcards.keys():
            if "sample" in wildcards.keys():
                return design[
                    (design["Sample_id"] == str(wildcards.sample)) & 
                    (design["Version"] == str(wildcards.version)) &
                    (design["Manip"] == str(wildcards.manip))
                ].Baseline.to_list()[0]
    
    # Same remark with manip.
    if "version" in wildcards.keys():
        if "sample" in wildcards.keys():
            return design[
                (design["Sample_id"] == str(wildcards.sample)) & 
                (design["Version"] == str(wildcards.version))
            ].Baseline.to_list()[0]

    # Same remark with panel version.
    if "sample" in wildcards.keys():
        return design[
            (design["Sample_id"] == str(wildcards.sample))
        ].Baseline.to_list()[0]

    raise ValueError("Missing wildcards values.")


# def get_wbc(wildcards):
#     return link_sample_baseline[wildcards.sample]["wbc"]

def get_wbc(wildcards):
    # Optimal case.
    if "nb" in wildcards.keys():
        if "manip" in wildcards.keys():
            if "version" in wildcards.keys():
                if "sample" in wildcards.keys():
                    return design[
                        (design["Sample_id"] == str(wildcards.sample)) & 
                        (design["Version"] == str(wildcards.version)) &
                        (design["Manip"] == str(wildcards.manip)) &
                        (design["NB"] == str(wildcards.nb))
                    ].WBC.to_list()[0]

    # There is only one baseline per replicate, so
    # missing "nb" should not be a problem.
    if "manip" in wildcards.keys():
        if "version" in wildcards.keys():
            if "sample" in wildcards.keys():
                return design[
                    (design["Sample_id"] == str(wildcards.sample)) & 
                    (design["Version"] == str(wildcards.version)) &
                    (design["Manip"] == str(wildcards.manip))
                ].WBC.to_list()[0]

    raise ValueError("Missing wildcards values: get_wbc requires 'manip'.")


# def get_ctc(wildcards):
#     return link_sample_baseline[wildcards.sample]["ctc"]

def get_ctc(wildcards):
    if "nb" in wildcards.keys():
        if "manip" in wildcards.keys():
            if "version" in wildcards.keys():
                if "sample" in wildcards.keys():
                    return design[
                        (design["Sample_id"] == str(wildcards.sample)) & 
                        (design["Version"] == str(wildcards.version)) &
                        (design["Manip"] == str(wildcards.manip)) &
                        (design["NB"] == str(wildcards.nb))
                    ].WBC.to_list()[0]
    
    raise ValueError("Missing wildcards values: get_ctc requires 'nb', 'manip', 'version', and 'sample'.")


# def get_trio(wildcards):
#     return {
#         "tumor": link_sample_baseline[wildcards.sample]["ctc"],
#         "wbc": link_sample_baseline[wildcards.sample]["wbc"],
#         "normal": link_sample_baseline[wildcards.sample]["baseline"],
#         "fasta": config["ref"]["fasta"],
#     }


def get_trio(wildcards):
    if "nb" in wildcards.keys():
        return {
            "tumor": get_ctc(wildcards),
            "normal": get_baseline(wildcards),
            "fasta": config["ref"]["fasta"],
        }
    return {
        "tumor": get_wbc(wildcards),
        "normal": get_baseline(wildcards),
        "fasta": config["ref"]["fasta"],
    }
    


# def get_trio_wbc(wildcards):
#     return {
#         "ctc": link_sample_baseline[wildcards.sample]["ctc"],
#         "tumor": link_sample_baseline[wildcards.sample]["wbc"],
#         "normal": link_sample_baseline[wildcards.sample]["baseline"],
#         "fasta": config["ref"]["fasta"],
#     }


# def get_hc(wildcards):        
#     return {
#         "bam": link_sample_baseline[wildcards.sample][wildcards.status],
#         "fasta": config["ref"]["fasta"],
#     }

def get_hc(wildcards):
    if "nb" in wildcards.keys():
        return {
            "bam": get_ctc(wildcards),
            "fasta": config["ref"]["fasta"],
        }
    if "manip" in wildcards.keys():
        return {
            "bam": get_wbc(wildcards),
            "fasta": config["ref"]["fasta"],
        }
    return {
        "bam": get_baseline(wildcards),
        "fasta": config["ref"]["fasta"],
    }


# def get_ensembl_vep_hc(wildcards):
#     subdir = None
#     if wildcards.status == "baseline":
#         subdir = "Baseline"
#     elif wildcards.status == "wbc":
#         subdir = "WBC"
#     elif wildcards.status == "ctc":
#         subdir = "HC_CTC"
#     else:
#         raise ValueError(f"Uknown status: {wildcards.status}")

#     return {
#         "cache": config["ref"]["vep"],
#         "fasta": "resources/GRCh38.fasta",
#         "vcf": f"data_output/{subdir}/{wildcards.sample}.vcf.gz",
#         "vcf_tbi": f"data_output/{subdir}/{wildcards.sample}.vcf.gz.tbi",
#     }

def get_ensembl_vep_hc(wildcards):
    if "nb" in wildcards.keys():
       return {
            "cache": config["ref"]["vep"],
            "fasta": "resources/GRCh38.fasta",
            "vcf": f"data_output/HC_CTC/{wildcards.sample}_{wildcards.version}_{wildcards.manip}_{wildcard.nb}.vcf.gz",
            "vcf_tbi": f"data_output/HC_CTC/{wildcards.sample}_{wildcards.version}_{wildcards.manip}_{wildcard.nb}.vcf.gz.tbi",
        }
    if "manip" in wildcards.keys():
        return {
            "cache": config["ref"]["vep"],
            "fasta": "resources/GRCh38.fasta",
            "vcf": f"data_output/HC_WBC/{wildcards.sample}_{wildcards.version}_{wildcards.manip}.vcf.gz",
            "vcf_tbi": f"data_output/HC_WBC/{wildcards.sample}_{wildcards.version}_{wildcards.manip}.vcf.gz.tbi",
        }
    return {
        "cache": config["ref"]["vep"],
        "fasta": "resources/GRCh38.fasta",
        "vcf": f"data_output/HC_Baseline/{wildcards.sample}.vcf.gz",
        "vcf_tbi": f"data_output/HC_Baseline/{wildcards.sample}.vcf.gz.tbi",
    }


# link_bams, samples_list, link_sample_baseline, baseline_sample_list, wbc_sample_list = parse_design(design.copy())

# sample_baseline_table = pandas.DataFrame.from_dict(link_sample_baseline, orient="index")
# sample_baseline_table.set_index(["baseline", "wbc"], inplace=True)
# logging.info(
#     f"First 20 lines of fastq correspondancies: \n{sample_baseline_table.head(20)}"
# )

# replicate_list = list(set(design["Replicate"]))
# version_list = list(set(design["Version"]))
# manip_list = list(set(design["Manip"]))
# status_list = list(set(design["Status"]))
# tmp = os.environ.get("BIGR_DEFAULT_TMP", "tmp")
# mutect_dir_list=["mutect2", "mutect2_wbc"]

ctc_list = list(set([
    f"{s}_{v}_{m}_{n}" for s, v, m, n in zip(
        design.Sample_id,
        design.Version,
        design.Manip,
        design.NB,
    )
]))

wbc_list = list(set(
    f"{s}_{v}_{m}" for s, v, m in zip(
        design.Sample_id,
        design.Version,
        design.Manip,
    )
))

baseline_list = list(set(design.Sample_id.to_list()))
status_list = ["CTC", "WBC", "Baseline"]
hc_list = [f"HC_{status}" for status in status_list]
mutect_list = [f"Mutect_{status}" for status in status_list[0:-1]]
tmp = os.environ.get("BIGR_DEFAULT_TMP", "tmp")

wildcard_constraints:
    sample=r"|".join(list(set(design.Sample_id.to_list()))),
    version=r"|".join(list(set(design.Version.to_list()))),
    manip=r"|".join(list(set(design.Manip.to_list()))),
    nb=r"|".join(list(set(design.NB.to_list()))),
    ctc_sample=r"|".join(ctc_list),
    wbc_sample=r"|".join(wbc_list),
    baseline_sample=r"|".join(baseline_list),
    # sample=r"|".join(samples_list + baseline_sample_list + wbc_sample_list),
    # status=r"|".join(status_list),
    # mutect_dir=r"|".join(mutect_dir_list),
    # sample_wbc=r"|".join(wbc_sample_list),
