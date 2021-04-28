"""
This file contains usefull functions for file managment:

* Search files (fastq, vcf)
* get an index from a basename (fasta, vcf, bcf, bam, sam, cram)
* get a dictionnary from a basename (fasta)
"""

import pandas

from pathlib import Path

from os.path import commonprefix
from typing import Any, Callable, Optional, Union

FilePathType = Union[str, Path]

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
    print("Searching in: {} for files ending with {}".format(dirpath, ext))
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


def search_vcf_files(dirpath: FilePathType) -> list[str]:
    return search_files(dirpath, ext=("vcf", "vcf.gz", "bcf"))

def search_vcf(dirpath: FilePathType) -> list[str]:
    return search_vcf_files(dirpath)

def search_fastq_files(dirpath: FilePathType) -> list[str]:
    return search_files(dirpath, ext=("fastq", "fq", "fastq.gz", "fq.gz"))

def search_fastq_pairs(dirpath: FilePathType) -> dict[str, dict[str, str]]:
    return {
        commonprefix([Path(r1).name, Path(r2).name]): {"Upstream_file": r1, "Downstream_file": r2}
        for r1, r2 in zip(*[iter(search_fastq_files(dirpath))]*2)
    }

def search_mapping(dirpath: FilePathType) -> list[str]:
    return search_files(dirpath, ext=("sam", "bam", "cram"))

def search_bam(dirpath: FilePathType) -> list[str]:
    return search_mapping(dirpath)

def search_fasta(dirpath: FilePathType) -> list[str]:
    return search_files(dirpath, ext=("fa", "fasta", "fa.gz", "fasta.gz"))

def search_fa(dirpath: FilePathType) -> list[str]:
    return search_fasta(dirpath)

def read_design(design_path: FilePathType) -> pandas.DataFrame:
    return pandas.read_csv(design_path, sep="\t", header=0)

def build_design(
        dirpath: FilePathType,
        search_fn: Callable[FilePathType, dict[str, dict[str, str]]]) \
        -> pandas.DataFrame:
    design = (pandas.DataFrame.from_dict(search_fn(dirpath), orient="index")
                              .reset_index()
                              .rename(columns={"index": "Sample_id"}))
    design["Sample_id"] = [x.strip(".") for x in design["Sample_id"]]
    design.to_csv("design.tsv", sep="\t", header=True, index=False)
    return design


def build_config(
        default_config: dict[str, Any],
        default_config_path: FilePathType) -> dict[str, Any]:
    if isinstance(default_config_path, str):
        default_config_path = Path(default_config_path)

    if default_config_path.exists():
        print("Configfile provided.")
    else:
        print("Missing config file, falling back to default parameters")
