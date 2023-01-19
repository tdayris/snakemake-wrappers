import datetime
import logging
import os
import pandas
import sys
import functools
from pathlib import Path

workflow_source_dir = Path(snakemake.workflow.srcdir(".."))
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


container: "docker://continuumio/miniconda3:4.4.10"


localrules:
    bigr_copy,


ruleorder: gatk_filter_mutect_calls > tabix_index
ruleorder: mutect2_somatic > tabix_index
ruleorder: gatk_variant_filtration > pbgzip_compress
ruleorder: gatk_variant_filtration > tabix_index
ruleorder: gleaves_compatibility > tabix_index
ruleorder: gleaves_compatibility > pbgzip_compress


########################
### Load environment ###
########################

default_config = read_yaml(workflow_source_dir / "config.hg38.yaml")


configfile: get_config(default_config)


wrapper_prefix = workflow_source_dir / ".." / ".."

design = get_design(os.getcwd(), search_fastq_somatic)
design.dropna(inplace=True)
design.index = design["Sample_id"]

########################
### Global variables ###
########################

tmp = os.getenv("BIGR_DEFAULT_TMP", "tmp")
sample_list = design["Sample_id"].tolist()
streams = ["1", "2"]
status_list = ["normal", "tumor"]
content = ["snp", "indel"]
cleaning_status = ["raw", "cleaned"]


fastq_links = link_fq_somatic(
    sample_names=sample_list,
    n1_paths=design.Upstream_file_normal,
    t1_paths=design.Upstream_file_tumor,
    n2_paths=design.Downstream_file_normal,
    t2_paths=design.Downstream_file_tumor,
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
