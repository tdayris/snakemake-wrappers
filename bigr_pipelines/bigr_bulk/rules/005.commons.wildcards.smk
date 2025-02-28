"""
This snakefile contains python functions for:
* Listing possible values taken by every single snakemake wildcard
* Constraining wildcards
"""

import itertools

sample_list = design["Sample_id"]
genome_builds = ["GRCh38", "GRCm39"]
genome_releases = ["108"]
stream_list = ["1", "2"]
sample_streams = [
    ".".join(sample_stream) 
    for sample_stream in itertools.product(sample_list, stream_list)
]

#############################
### RNA-Seq related lists ###
#############################
rna_sample_list = design[design["Protocol"].str.lower.isin(("rnaseq", "unknown"))]["Sample_id"]


wildcard_constraints:
    sample=r"|".join(sample_list),
    genome_build=r"|".join(genome_builds),
    genome_releases=r"|".join(genome_releases),
    stream=r"|".join(stream_list),
    sample_stream=r"|".join(sample_streams),