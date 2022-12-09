#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
This file contains usefull functions for file managment:

* Search files (fastq, vcf)
* get an index from a basename (fasta, vcf, bcf, bam, sam, cram)
* get a dictionnary from a basename (fasta)
"""

import logging
import pandas

from pathlib import Path

from os.path import basename
from typing import Any, Callable, Dict, List, Optional, Union

from write_yaml import *

FilePathType = Union[str, Path]


def remove_suffixes(word: str, suffixes: List[str]) -> str:
    """
    From a given word (usually a path) remove the provided suffixe
    (usually an extension)
    """
    for suffix in suffixes:
        if word.endswith(suffix):
            return word[:-len(suffix)]
    return word


def get_fasta_index_from_genome_path(genome_path: str) -> str:
    """
    From a fasta path, return it's fai without checking for existence (faster)
    """
    return genome_path + ".fai"


def get_fai(genome_path: str) -> str:
    """
    From a fasta path, return it's fai without checking for existence (faster)
    """
    return get_fasta_index_from_genome_path(genome_path)

def get_fasta_dict_from_genome_path(genome_path: str) -> str:
    """
    From a fasta path, return it's dictionnary without checking for
    existence (faster)
    """
    return ".".join(genome_path.split(".")[:-1]) + ".dict"


def get_dict(genome_path: str) -> str:
    """
    From a fasta path, return it's dictionnary without checking for
    existence (faster)
    """
    return get_fasta_dict_from_genome_path(genome_path)


def get_bam_index_from_bam_path(bam_path: str) -> str:
    """
    From a bam/sam path, return the BAI index
    """
    return f"{bam_path}.bai"


def get_bai(bam_path: str) -> str:
    """
    From a bam/sam path, return the BAI index
    """
    return get_bam_index_from_bam_path(bam_path)


def get_vcf_tbi_from_vcf_path(vcf_path: str) -> str:
    """
    From a vcf path, return it's tbi index without checking for
    existence (faster)
    """
    return vcf_path + ".tbi"


def get_tbi(vcf_path: str) -> str:
    """
    From a vcf path, return it's tbi index without checking for
    existence (faster)
    """
    return get_vcf_tbi_from_vcf_path(vcf_path)


def search_files(dirpath: FilePathType, ext: Optional[str] = None) -> List[str]:
    """
    Within a given directory, search all files. If an extension is provided,
    return only the ones with matching suffixe. If a path leads to a directory,
    this function is called recursively.
    """
    logging.info("Searching in: {} for files ending with {}".format(dirpath, ext))
    if isinstance(dirpath, str):
        dirpath = Path(dirpath)

    # Case user provided a path to a file instead of a path to a directory
    if not dirpath.is_dir():
        dirpath = dirpath.parent

    for filepath in dirpath.iterdir():
        if filepath.is_dir() and filepath.name != ".snakemake":
            # Recursive search
            yield from search_files(filepath, ext)
        elif (ext is None) or filepath.name.endswith(ext):
            yield str(filepath.absolute())


def search_vcf_files(dirpath: FilePathType) -> Dict[str, str]:
    """
    Within a given directory, search all vcf files. If a path leads to a
    directory, the subfunction is called recursively.
    """
    suffixes = ("vcf", "vcf.gz", "bcf")
    return {
        remove_suffixes(basename(vcf), suffixes): {"Upstream_file": vcf}
        for vcf in search_files(dirpath, ext=suffixes)
    }


def search_vcf(dirpath: FilePathType) -> Dict[str, str]:
    """
    Within a given directory, search all vcf files. If a path leads to a
    directory, the subfunction is called recursively.
    """
    return search_vcf_files(dirpath)


def search_fastq_files(dirpath: FilePathType) -> Dict[str, str]:
    """
    Within a given directory, search all fastq files. If a path leads to a
    directory, the subfunction is called recursively.
    """
    suffixes = (".fastq", ".fq", ".fastq.gz", ".fq.gz")
    return {
        remove_suffixes(basename(fastq), suffixes): fastq
        for fastq in sorted(search_files(dirpath, ext=suffixes))
    }


def search_fastq_files_dict(dirpath: FilePathType) -> Dict[str, Dict[str, str]]:
    """
    Within a given directory, search all fastq files. If a path leads to a
    directory, the subfunction is called recursively.
    """
    suffixes = (".fastq", ".fq", ".fastq.gz", ".fq.gz")
    return {
        remove_suffixes(basename(fastq), suffixes): {"Upstream_file": fastq}
        for fastq in sorted(search_files(dirpath, ext=suffixes))
    }


def search_fastq_somatic(dirpath: FilePathType) -> Dict[str, str]:
    """
    Within a given directory, search all vcf files. If a path leads to a
    directory, the subfunction is called recursively.

    It returns fastq files four by four, according to alphanumerical order.
    """
    suffixes = ("fastq", "fq", "fastq.gz", "fq.gz")
    return {
        remove_suffixes(basename(t1), suffixes): {
            "Upstream_file_tumor": t1,
            "Downstream_file_tumor": t2,
            "Upstream_file_normal": n1,
            "Downstream_file_normal": n2
        } for n1, n2, t1, t2 in zip(*[iter(search_fastq_files(dirpath).values())]*4)
    }


def search_fastq_trio(dirpath: FilePathType) -> Dict[str, str]:
    """
    Within a given directory, search all vcf files. If a path leads to a
    directory, the subfunction is called recursively.

    It returns fastq files six by six, according to alphanumerical order.
    """
    suffixes = ("fastq", "fq", "fastq.gz", "fq.gz")
    fqiter = zip(*[iter(search_fastq_files(dirpath))]*6)
    return {
        remove_suffixes(basename(n1), suffixes): {
            "Tumor_upstream_file": n1,
            "Tumor_downstream_file": n2,
            "Father_upstream_file": f1,
            "Father_downstream_file": f2,
            "Mother_upstream_file": m1,
            "Mother_downstream_file": m2,
        } for n1, n2, f1, f2, m1, m2 in fqiter
    }


def search_fastq_pairs(dirpath: FilePathType) -> Dict[str, Dict[str, str]]:
    """
    Within a given directory, search all vcf files. If a path leads to a
    directory, the subfunction is called recursively.

    It returns fastq files as pairs, according to alphanumerical order.
    """
    suffixes = (".fastq", ".fq", ".fastq.gz", ".fq.gz")
    return {
        remove_suffixes(basename(r1), suffixes): {
            "Upstream_file": r1, "Downstream_file": r2
        } for r1, r2 in zip(*[iter(sorted(search_fastq_files(dirpath).values()))]*2)
    }


def search_mapping(dirpath: FilePathType) -> Dict[str, str]:
    """
    Within a given directory, search all bam files. If a path leads to a
    directory, the subfunction is called recursively.
    """
    suffixes = ("sam", "bam", "cram")
    return {
        remove_suffixes(basename(mapping), suffixes): mapping
        for mapping in search_files(dirpath, ext=suffixes)
    }


def search_bam(dirpath: FilePathType) -> Dict[str, str]:
    """
    Within a given directory, search all bam files. If a path leads to a
    directory, the subfunction is called recursively.
    """
    return search_mapping(dirpath)


def search_fasta(dirpath: FilePathType) -> List[str]:
    """
    Within a given directory, search all fasta files. If a path leads to a
    directory, the subfunction is called recursively.
    """
    return search_files(dirpath, ext=("fa", "fasta", "fa.gz", "fasta.gz"))


def search_fa(dirpath: FilePathType) -> List[str]:
    """
    Within a given directory, search all fasta files. If a path leads to a
    directory, the subfunction is called recursively.
    """
    return search_fasta(dirpath)


def read_design(design_path: FilePathType) -> pandas.DataFrame:
    """
    Load provided design file as a pandas DataFrame.
    """
    return pandas.read_csv(design_path, sep="\t", header=0)


def build_design(
        dirpath: FilePathType,
        search_func: Callable[FilePathType, Dict[str, Dict[str, str]]]) \
        -> pandas.DataFrame:
    """
    Search for files with the provided function within the given directory,
    then return a design as a pandas DataFrame.
    """
    design = (pandas.DataFrame.from_dict(search_func(dirpath), orient="index")
                              .reset_index()
                              .rename(columns={"index": "Sample_id"}))
    design["Sample_id"] = [x.strip(".") for x in design["Sample_id"]]
    design.to_csv("design.tsv", sep="\t", header=True, index=False)
    return design


def get_design(
        dirpath: FilePathType,
        search_func: Callable[FilePathType, Dict[str, Dict[str, str]]]) \
        -> pandas.DataFrame:
    """
    From a given path, search for a file named 'design.tsv'. If this file is
    missing, then this function searches for files with the given extension
    within the provided directory, then builds, saves and returns a pandas
    DataFrame corresponding to that design.
    """
    design = None
    if (design_path := Path("design.tsv")).exists():
        logging.info("Design file available: %s", str(design_path.absolute()))
        design = read_design(design_path)
    else:
        logging.info("Design file not found, looking for avilable input files")
        design = build_design(dirpath, search_func)
    return design


def get_config(default_config: Dict[str, Any]) -> str:
    """
    From a given path, seach for a function named "config.yaml". If this file
    is missing, this function saves a copy of the default configuration.
    """
    if (config_path := Path("config.yaml")).exists():
        logging.info("Config file available: %s", str(config_path.absolute()))
    else:
        logging.info(
            "Generic config file not found, falling back to default arguments"
        )
        if isinstance(default_config, str):
            write_yaml(output_yaml=config_path, data=read_yaml(default_config))
        else:
            write_yaml(output_yaml=config_path, data=default_config)
    return str(config_path.absolute())


def design_config(
        default_config: Dict[str, Any],
        dirpath: FilePathType,
        search_func: Callable[FilePathType, Dict[str, Dict[str, str]]]) \
        -> List[Union[pandas.DataFrame, str]]:
    """
    Shortcut to build/load both design and config at once.
    """
    return [get_design(dirpath, search_func), get_config(default_config)]
