import datetime
import logging
import os
import pandas
import sys
import functools
from pathlib import Path

# workflow_source_dir = Path(snakemake.workflow.srcdir(".."))
common = str(workflow_source_dir / ".." / "common" / "python")
sys.path.append(common)

from file_manager import *
from files_linker import *
from write_yaml import *
from messages import *
from reservation import *
from snakemake.utils import min_version

min_version("6.0")

logging.basicConfig(
    filename="snakemake.variant_calling_somatic.log", filemode="w", level=logging.DEBUG
)

default_container_path = Path(
    workflow_source_dir / ".." / ".." / ".." / "singularity" / "mambaforge_4.14.0-0.sif"
)
container_path = (
    str(default_container_path)
    if default_container_path.exists
    else "docker://continuumio/miniconda3:4.4.10"
)


container: container_path


localrules:
    bigr_copy,


ruleorder: gatk_filter_mutect_calls > tabix_index
ruleorder: mutect2_somatic > tabix_index
ruleorder: gatk_variant_filtration > pbgzip_compress
ruleorder: gatk_variant_filtration > tabix_index
ruleorder: gleaves_compatibility > tabix_index
ruleorder: gleaves_compatibility > pbgzip_compress
ruleorder: gatk_hard_filtering > tabix_index
ruleorder: gatk_apply_baserecalibrator > sambamba_index_bam


########################
### Load environment ###
########################

default_config = read_yaml(workflow_source_dir / "config.hg38.yaml")


configfile: get_config(default_config)


wrapper_prefix = workflow_source_dir / ".. " / ".."

design = get_design(os.getcwd(), search_fastq_somatic)
design.dropna(inplace=True)
design.index = design["Sample_id"]

########################
### Global variables ###
########################

tmp = os.getenv("BIGR_DEFAULT_TMP", "tmp")
sample_list = design["Sample_id"].tolist()
streams = ["1", "2"]
status_list = (
    ["normal", "tumor"]
    if "Upstream_file_normal" in design.columns.tolist()
    else ["tumor"]
)
content = ["snp", "indel"]
cleaning_status = ["raw", "cleaned"]


fastq_links = link_fq_somatic(
    sample_names=sample_list,
    n1_paths=design.get("Upstream_file_normal"),
    t1_paths=design.get("Upstream_file_tumor"),
    n2_paths=design.get("Downstream_file_normal"),
    t2_paths=design.get("Downstream_file_tumor"),
    prefix="data_input",
)

# Handle mouse missing databases
last_vcf = (
    "bigr/cancer_gene_census/{sample}.vcf"
    if config["reference"]["ncbi_build"] != "mm10"
    else "snpsift/dbsnp/{sample}.vcf"
)


wildcard_constraints:
    sample=r"|".join(sample_list),
    stream=r"|".join(streams),
    status=r"|".join(status_list),
    content=r"|".join(content),
    cleaning=r"|".join(cleaning_status),


#########################
### Herlper functions ###
#########################


def get_mutect2_input(wildcards) -> Dict[str, str]:
    """
    Given somatic/germline design file, return the correct
    tumor bam + optional normal bam pairs
    """
    base = {
        "fasta": config["reference"]["fasta"],
        "fasta_idx": config["reference"]["fasta_index"],
        "fasta_dict": config["reference"]["fasta_dict"],
        "germline": config["reference"]["af_only"],
        "germline_tbi": config["reference"]["af_only_tbi"],
        "intervals": config["reference"]["capture_kit_bed"],
    }
    if "Upstream_file_normal" in design.columns.tolist():
        base["map"] = f"sambamba/markdup/{wildcards.sample}_normal.bam"
        base["map_idx"] = f"sambamba/markdup/{wildcards.sample}_normal.bam.bai"
        base["tumor"] = f"sambamba/markdup/{wildcards.sample}_tumor.bam"
        base["tumor_idx"] = f"sambamba/markdup/{wildcards.sample}_tumor.bam.bai"
    else:
        base["map"] = f"sambamba/markdup/{wildcards.sample}_tumor.bam"
        base["map_idx"] = f"sambamba/markdup/{wildcards.sample}_tumor.bam.bai"

    return base


def get_mutect2_args(wildcards) -> str:
    """
    Return Mutect2 optional arguments with sample decoration
    in case of somatic variant calling
    """
    base = config["gatk"].get("mutect2", "")
    if "Upstream_file_normal" in design.columns.tolist():
        base += f"--tumor-sample {wildcards.sample}_tumor "
        base += f"--normal-sample {wildcards.sample}_normal "
    return base


def get_filter_mutect2_input(wildcards) -> Dict[str, str]:
    """
    No tumor contamination estimates is possible without
    normal calls.
    """
    base = {
        "vcf": "mutect2/call/{sample}.vcf.gz",
        "vcf_tbi": get_tbi("mutect2/call/{sample}.vcf.gz"),
        "ref": config["reference"]["fasta"],
        "ref_index": config["reference"]["fasta_index"],
        "ref_dict": config["reference"]["fasta_dict"],
        "bam": "sambamba/markdup/{sample}_tumor.bam",
        "bam_index": get_bai("sambamba/markdup/{sample}_tumor.bam"),
        "f1r2": "gatk/orientation_model/{sample}/{sample}.artifacts-prior.tar.gz",
    }

    if "Upstream_file_normal" in design.columns.tolist():
        base["contamination"] = "summary/{sample}_calculate_contamination.table"

    return base


def targets(wildcards):
    base = {
        "bam": expand(
            "data_output/BAM/{sample}_{status}.bam",
            sample=sample_list,
            status=status_list,
        ),
        "bam_md5": expand(
            "data_output/BAM/{sample}_{status}.bam.md5",
            sample=sample_list,
            status=status_list,
        ),
        "bai": expand(
            "data_output/BAM/{sample}_{status}.bam.bai",
            sample=sample_list,
            status=status_list,
        ),
        "mapping_QC": "data_output/MultiQC/MappingQC.html",
    }

    if config.get("steps", {}).get("calling", True):
        base["vcf"] = expand(
            "data_output/VCF/{sample}.vcf.gz",
            sample=sample_list,
        )
        base["vcf_tbi"] = expand(
            "data_output/VCF/{sample}.vcf.gz.tbi",
            sample=sample_list,
        )
        base["vcf_md5"] = expand(
            "data_output/VCF/{sample}.vcf.gz.md5",
            sample=sample_list,
        )
        base["tsv"] = expand(
            "data_output/TSV/{sample}.tsv",
            sample=sample_list,
        )
        base["xslx"] = expand(
            "data_output/TSV/{sample}.tsv",
            sample=sample_list,
        )
        base["calling_qc"] = "data_output/MultiQC/Somatic_Variant_Calling.html"

    if config.get("steps", {}).get("cnv", False):
        if "Upstream_file_normal" in design.columns.tolist():
            base["cnv"] = expand(
                "data_output/CNV/{sample}.tsv",
                sample=sample_list,
            )
        else:
            logging.warning(
                "CNV are not analyzed without Normal/Tumor somatic calling."
            )

    if config.get("steps", {}).get("tmb", False):
        base["tmb"] = "data_output/TMB.tsv"

    if config.get("steps", {}).get("msi", False):
        if "Upstream_file_normal" in design.columns.tolist():
            base["msi"] = "data_output/MSI.tsv"
            base["calling_qc"] = "data_output/MultiQC/Somatic_Variant_Calling.html"
        else:
            logging.warning(
                "MSI cannot be analyzed without Normal/Tumor somatic calling."
            )

    return base
