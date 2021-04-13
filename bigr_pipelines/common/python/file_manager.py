"""
This file contains usefull functions for file managment:

* Search files (fastq, vcf)
* get an index from a basename (fasta, vcf, bcf, bam, sam, cram)
* get a dictionnary from a basename (fasta)
"""

from pathlib import Path

from os.path import commonprefix
from typing import Optional, Union

FilePathType = Union[str, Path]

def get_fasta_index_from_genome_path(genome_path: str) -> str:
    return genome_path + ".fai"

def get_fasta_dict_from_genome_path(genome_path: str) -> str:
    return ".".join(genome_path.split(".")[:-1]) + ".dict"


def get_vcf_tbi_from_vcf_path(vcf_path: str) -> str:
    return db_path + ".tbi"

def search_files(dirpath: FilePathType, ext: Optional[str] = None) -> list[str]:
    result = []
    if isinstance(dirpath, str):
        dirpath = Path(dirpath)

    for filepath in dirpath.iterdir():
        if ext is None:
            result.append(filepath)
        elif filepath.name.endswith(ext):
            result.append(filepath)


def search_vcf_files(dirpath: FilePathType) -> list[str]:
    return search_files(dirpath, ext=("vcf", "vcf.gz", "bcf"))

def search_fastq_files(dirpath: FilePathType) -> list[str]:
    return search_files(dirpath, ext=("fastq", "fq", "fastq.gz", "fq.gz"))

def search_fastq_pairs(dirpath: FilePathType) -> dict[str, list[str]]:
    return {
        commonprefix([Path(r1).name, Path(r2).name]): [r1, r2]
        for r1, r2 in zip(*[iter(search_fastq_files(dirpath))]*2)
    }

def search_mapping(dirpath: FilePathType) -> list[str]:
    return search_files(dirpath, ext=("sam", "bam", "cram"))

def search_fasta(dirpath: FilePathType) -> list[str]:
    return search_files(dirpath, ext=("fa", "fasta", "fa.gz", "fasta.gz"))

def read_design(design_path: FilePathType) -> pandas.DataFrame:
    return pandas.read_csv(design_path, sep="\t", header=0)
