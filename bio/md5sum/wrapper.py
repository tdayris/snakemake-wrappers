#!/usr/bin/python3.8
# -*- coding: utf-8 -*-

"""Snakemake wrapper for md5sum"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2020, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

from snakemake.shell import shell
log = snakemake.log_fmt_shell(stdout=False, stderr=True)


def read_md5(file_path: str) -> str:
    """
    Read a MD5sum output file and return the hash value
    """
    with open(file_path) as hashfile:
        hash_value, *_ = next(hashfile).split(" ")

    return hash_value


def compare_md5(left: str, right: str) -> None:
    """
    Compare two md5 signatures
    """
    if not left == right:
        raise ValueError(
            "MD5 hash was not the expected one: "
            "got     : {}\nexpected: {}".format(left, right)
        )


shell("md5sum {snakemake.input.file} > {snakemake.output.hash} {log}")


if (hash_value := snakemake.params.get("hash_value")) is not None:
    compare_md5(hash_value, read_md5(snakemake.output[0]))
elif (hash_file := snakemake.input.get("hash")) is not None:
    compare_md5(read_md5(hash_file), read_md5(snakemake.output[0]))
