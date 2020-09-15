#!/usr/bin/python3.7
# -*- coding: utf-8 -*-

"""Snakemake wrapper md5sum checking"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2019, Thibault Dayris"
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


expected = ""
if snakemake.params.get("hash_value", None) is not None:
    expected = snakemake.params["hash_value"]
elif "hash" in snakemake.input.keys():
    expected = read_md5(snakemake.input["hash"])
else:
    raise ValueError(
        "Could not determine expected hash value for MD5 comparison"
    )

shell("md5sum {snakemake.input.file} > {snakemake.output} {log}")

obtained = read_md5(snakemake.output[0])

if not expected == obtained:
    raise ValueError(
        "MD5 hash was not the expected one: "
        "got {}, expected {}".format(obtained, expected)
    )
