import datetime
import logging
import os
import pandas
import sys
from pathlib import Path

worflow_source_dir = Path(snakemake.workflow.srcdir("."))
common = str(worflow_source_dir / ".." / "common" / "python")
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


default_config = read_yaml(worflow_source_dir / "config.hg38.yaml")


configfile: get_config(default_config)


design = get_design(os.getcwd(), search_fastq_somatic)
# design = design.head(2).tail(1)
design.dropna(inplace=True)
# print(design)

design.index = design["Sample_id"]
# design.drop(index="s070", inplace=True)


wildcard_constraints:
    sample=r"|".join(design["Sample_id"]),
    stream=r"1|2|R1|R2",
    status=r"normal|tumor",
    content=r"snp|indel",


fastq_links = link_fq_somatic(
    sample_names=design.Sample_id,
    n1_paths=design.Upstream_file_normal,
    t1_paths=design.Upstream_file_tumor,
    n2_paths=design.Downstream_file_normal,
    t2_paths=design.Downstream_file_tumor,
)


ruleorder: fix_annotation_for_gatk > pbgzip_compress
ruleorder: gatk_variant_filtration > pbgzip_compress
