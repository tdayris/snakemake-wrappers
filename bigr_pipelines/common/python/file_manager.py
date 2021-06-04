"""
This file contains usefull functions for file managment:

* Search files (fastq, vcf)
* get an index from a basename (fasta, vcf, bcf, bam, sam, cram)
* get a dictionnary from a basename (fasta)
"""

import logging
import pandas
import os
import shutil

from pathlib import Path

from os.path import commonprefix, basename
from typing import Any, Callable, Optional, Union

from write_yaml import *

FilePathType = Union[str, Path]


def remove_suffixes(word: str, suffixes: list[str]) -> str:
    for suffix in suffixes:
        if word.endswith(suffix):
            return word[:-len(suffix)]
    return word


def get_fasta_index_from_genome_path(genome_path: str) -> str:
    return genome_path + ".fai"


def get_fai(genome_path: str) -> str:

    return get_fasta_index_from_genome_path(genome_path)

def get_fasta_dict_from_genome_path(genome_path: str) -> str:
    return ".".join(genome_path.split(".")[:-1]) + ".dict"


def get_dict(genome_path: str) -> str:
    return get_fasta_dict_from_genome_path(genome_path)


def get_vcf_tbi_from_vcf_path(vcf_path: str) -> str:
    return vcf_path + ".tbi"


def get_tbi(vcf_path: str) -> str:
    return get_vcf_tbi_from_vcf_path(vcf_path)


def search_files(dirpath: FilePathType, ext: Optional[str] = None) -> list[str]:
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


def search_vcf_files(dirpath: FilePathType) -> dict[str, str]:
    suffixes = ("vcf", "vcf.gz", "bcf")
    return {
        remove_suffixes(basename(vcf), suffixes): vcf
        for vcf in search_files(dirpath, ext=suffixes)
    }


def search_vcf(dirpath: FilePathType) -> dict[str, str]:
    return search_vcf_files(dirpath)


def search_fastq_files(dirpath: FilePathType) -> dict[str, str]:
    suffixes = (".fastq", ".fq", ".fastq.gz", ".fq.gz")
    return {
        remove_suffixes(basename(fastq), suffixes): fastq
        for fastq in sorted(search_files(dirpath, ext=suffixes))
    }


def search_fastq_somatic(dirpath: FilePathType) -> dict[str, str]:
    suffixes = ["fastq", "fq", "fastq.gz", "fq.gz"]
    return {
        remove_suffixes(basename(t1), suffixes): {
            "Tumor_upstream_file": t1,
            "Tumor_downstream_file": t2,
            "Normal_upstream_file": n1,
            "Normal_downstream_file": n2
        } for n1, n2, t1, t2 in zip(*[iter(search_fastq_files(dirpath))]*4)
    }


def search_fastq_pairs(dirpath: FilePathType) -> dict[str, dict[str, str]]:
    suffixes = (".fastq", ".fq", ".fastq.gz", ".fq.gz")
    return {
        remove_suffixes(basename(r1), suffixes): {
            "Upstream_file": r1, "Downstream_file": r2
        } for r1, r2 in zip(*[iter(sorted(search_fastq_files(dirpath).values()))]*2)
    }


def search_mapping(dirpath: FilePathType) -> dict[str, str]:
    suffixes = ("sam", "bam", "cram")
    return {
        remove_suffixes(basename(mapping), suffixes): mapping
        for mapping in search_files(dirpath, ext=suffixes)
    }


def search_bam(dirpath: FilePathType) -> dict[str, str]:
    return search_mapping(dirpath)


def search_fasta(dirpath: FilePathType) -> list[str]:
    return search_files(dirpath, ext=("fa", "fasta", "fa.gz", "fasta.gz"))


def search_fa(dirpath: FilePathType) -> list[str]:
    return search_fasta(dirpath)


def read_design(design_path: FilePathType) -> pandas.DataFrame:
    return pandas.read_csv(design_path, sep="\t", header=0)


def build_design(
        dirpath: FilePathType,
        search_func: Callable[FilePathType, dict[str, dict[str, str]]]) \
        -> pandas.DataFrame:
    design = (pandas.DataFrame.from_dict(search_func(dirpath), orient="index")
                              .reset_index()
                              .rename(columns={"index": "Sample_id"}))
    design["Sample_id"] = [x.strip(".") for x in design["Sample_id"]]
    design.to_csv("design.tsv", sep="\t", header=True, index=False)
    return design


def get_design(
        dirpath: FilePathType,
        search_func: Callable[FilePathType, dict[str, dict[str, str]]]) \
        -> pandas.DataFrame:
    design = None
    if (design_path := Path("design.tsv")).exists():
        logging.info("Design file available: %s", str(design_path.absolute()))
        design = read_design(design_path)
    else:
        logging.info("Design file not found, looking for avilable input files")
        design = build_design(dirpath, search_func)
    return design


def get_config(default_config: dict[str, Any]) -> str:
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
        default_config: dict[str, Any],
        dirpath: FilePathType,
        search_func: Callable[FilePathType, dict[str, dict[str, str]]]) \
        -> list[Union[pandas.DataFrame, str]]:
    return [get_design(dirpath, search_func), get_config(default_config)]
