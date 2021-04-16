"""
This file contains function to perform files mapping between original name
(hashed iRODS, biopath, ...) to normalized paths.
"""

from typing import Optional

def link_fq(
        sample_names: list[str],
        r1_paths: list[str],
        r2_paths: Optional[list[str]] = None
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
    if r2_paths is None:
        return {
            r1: f"reads/{sample}.fq.gz"
            for sample, r1 in zip(sample_names, r1_paths)
        }

    link_dict = {}
    for sample, r1, r2 in zip(sample_names, r1_paths, r2_paths):
        link_dict[f"reads/{sample}.1.fq.gz"] = r1
        link_dict[f"reads/{sample}.2.fq.gz"] = r2
    return link_dict
