"""
This file contains function to perform files mapping between original name
(hashed iRODS, biopath, ...) to normalized paths.
"""

import pandas

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
            f"reads/{sample}.fq.gz": r1
            for sample, r1 in zip(sample_names, r1_paths)
        }

    link_dict = {}
    for sample, r1, r2 in zip(sample_names, r1_paths, r2_paths):
        link_dict[f"reads/{sample}.1.fq.gz"] = r1
        link_dict[f"reads/{sample}.2.fq.gz"] = r2
    return link_dict


def link_fq_somatic(
        sample_names: list[str],
        n1_paths: list[str],
        t1_paths: list[str],
        n2_paths: Optional[list[str]] = None,
        t2_paths: Optional[list[str]] = None
    ) -> dict[str, dict[str, str]]:
    """
    Case r2 are provided:
    Build a dictionnary containing the following pairs:
    normal:
        original_r1_name: reads/normal/{sample}.1.fq.gz
        original_r2_name: reads/normal/{sample}.2.fq.gz
    tumor:
        original_r1_name: reads/tumor/{sample}.1.fq.gz
        original_r2_name: reads/tumor/{sample}.2.fq.gz

    Otherwise:
    Build a dictionnary containing the following fastq:
    tumor:
        original_name: reads/tumor/{sample}.fq.gz
    normal:
        original_name: reads/normal/{sample}.fq.gz
    """
    link_dict = {
        "normal": {},
        "tumor": {}
    }
    if n2_paths is None:
        for sample, n1 in zip(sample_names, n1_paths):
            link_dict["normal"][f"reads/normal/{sample}.fq.gz"] = n1
    else:
        for sample, n1, n2 in zip(sample_names, n1_paths, n2_paths):
            link_dict["normal"][f"reads/normal/{sample}.1.fq.gz"] = n1
            link_dict["normal"][f"reads/normal/{sample}.2.fq.gz"] = n2

    if t2_paths is None:
        for sample, t1 in zip(sample_names, t1_paths):
            result["tumor"][f"reads/tumor/{sample}.fq.gz"] = t1
    else:
        for sample, t1, t2 in zip(sample_names, t1_paths, t2_paths):
            link_dict["tumor"][f"reads/tumor/{sample}.1.fq.gz"] = t1
            link_dict["tumor"][f"reads/tumor/{sample}.2.fq.gz"] = t2

    return link_dict


def link_vcf(sample_names: list[str], files: list[str]) -> dict[str, str]:
    """
    Build a dictionnary linking a vcf path to its axpected path
    """
    return {
        f"data_input/calls/{sample}.vcf.gz": file
        for file, sample in zip(files, sample_names)
    }


def link_bed(design: pandas.DataFrame, bed: Optional[str] = None):
    """
    Build a dictionnary linking a sample name to a capturekit bed.

    Since more than one capturekig may be used in a sequencing run,
    this must be either implemented, or multiple instances of this
    pipelines are being run simultaneously: once per capture kit used.

    This is useless in a single-project single-captured analysis. So, yes,
    useless for data analysts, but usefull for piREST and automated pipelines
    """
    bed_to_sample_link = {}
    dataset_to_bed_path = {}
    if "datasetIdBED" in design.columns:
        for sample, sample_bed in zip(design["Sample_id", "datasetIdBED"]):
            bed_to_sample_link[sample] = dataset_to_bed_path.get(
                sample_bed,
                sample_bed if sample_bed is not None else bed
            )
    else:
        bed_to_sample_link = {
            sample: bed for sample in design["Sample_id"]
        }
    return bed_to_sample_link
