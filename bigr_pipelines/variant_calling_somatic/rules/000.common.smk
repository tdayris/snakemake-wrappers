import datetime
import logging
import os
import pandas
import sys
from pathlib import Path

workflow_source_dir = Path(snakemake.workflow.srcdir(".."))
common = str(workflow_source_dir / ".." / "common" / "python")
sys.path.append(common)

from file_manager import *
from files_linker import *
from write_yaml import *
from messages import *
from snakemake.utils import min_version

min_version("6.0")

logging.basicConfig(
    filename="snakemake.variant_calling_somatic.log", filemode="w", level=logging.DEBUG
)


container: "docker://continuumio/miniconda3:4.4.10"


localrules:
    bigr_copy,


ruleorder: sambamba_index_bam > sambamba_index
ruleorder: gatk_filter_mutect_calls > tabix_index
ruleorder: mutect2_somatic > tabix_index
ruleorder: create_snpeff_snpfit_data_input_dir > tabix_index
ruleorder: fix_annotation_for_gatk > pbgzip_compress
ruleorder: gatk_variant_filtration > pbgzip_compress


########################
### Load environment ###
########################

default_config = read_yaml(workflow_source_dir / "config.hg38.yaml")


configfile: get_config(default_config)


design = get_design(os.getcwd(), search_fastq_somatic)
design.dropna(inplace=True)
design.index = design["Sample_id"]

########################
### Global variables ###
########################

sample_list = design["Sample_id"].tolist()
streams = ["1", "2"]
statue = ["normal", "tumor"]
content = ["snp", "indel"]


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
    status=r"|".join(status),
    content=r"|".join(content),


############################
### Resource reservation ###
############################

# Memory and time reservation
def get_resources_per_gb(
    wildcards, input, attempt, multiplier: int = 1, base: int = 0
) -> int:
    """
    Return the amount of resources needed per GB of input.

    Parameters:
    * wildcards: Snakemake signature requires this parameter.
    * input: The input file list
    * attempt: The # of times the calling rule has been restarted
    * multiplier: An arbitrary multiplier

    Return:
    (int) The amount of resources needed (mb, minutes, etc)
    """
    return max(
        # Case there is 1gb or more in input
        # ((input.size // 10_000_000_000) * attempt * multiplier) + base,
        1,
        # Case there is less than 1gb in input
        (multiplier * attempt) + base,
    )


logging.info("Preparing memory calls...")
# Explicit time reservations
get_15min_per_attempt = functools.partial(get_resources_per_gb, multiplier=15)
get_45min_per_attempt = functools.partial(get_resources_per_gb, multiplier=45)
get_75min_per_attempt = functools.partial(get_resources_per_gb, multiplier=75)
get_20min_per_attempt = functools.partial(get_resources_per_gb, multiplier=20)
get_1h_per_attempt = functools.partial(get_resources_per_gb, multiplier=60)
get_2h_per_attempt = functools.partial(get_resources_per_gb, multiplier=60 * 2)
get_2h_per_attempt = functools.partial(get_resources_per_gb, multiplier=60 * 3)
get_5h_per_attempt = functools.partial(get_resources_per_gb, multiplier=60 * 5)
get_90min_per_attempt = functools.partial(get_resources_per_gb, multiplier=90)

# Explicit memory reservation
get_1gb_per_attempt = functools.partial(get_resources_per_gb, multiplier=1024)
get_2gb_per_attempt = functools.partial(get_resources_per_gb, multiplier=1024 * 2)
get_4gb_per_attempt = functools.partial(get_resources_per_gb, multiplier=1024 * 4)
get_6gb_per_attempt = functools.partial(get_resources_per_gb, multiplier=1024 * 6)
get_8gb_per_attempt = functools.partial(get_resources_per_gb, multiplier=1024 * 8)
get_10gb_per_attempt = functools.partial(get_resources_per_gb, multiplier=1024 * 10)
get_15gb_per_attempt = functools.partial(get_resources_per_gb, multiplier=1024 * 15)
get_20gb_per_attempt = functools.partial(get_resources_per_gb, multiplier=1024 * 20)
get_75gb_and_5gb_per_attempt = functools.partial(
    get_resources_per_gb, multiplier=1024 * 5, base=1024 * 75
)
get_20gb_and_10gb_per_attempt = functools.partial(
    get_resources_per_gb, multiplier=1024 * 10, base=1024 * 20
)
get_30gb_and_10gb_per_attempt = functools.partial(
    get_resources_per_gb, multiplier=1024 * 10, base=1024 * 30
)
