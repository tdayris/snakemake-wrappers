#!/usr/bin/python3.7
# -*- coding: utf-8 -*-

"""Snakemake wrapper for bash copy within the IGR's Flamingo cluster"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2019, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

import os.path as op

from snakemake.shell import shell
from tempfile import TemporaryDirectory


# Prepare logging
log = snakemake.log_fmt_shell(stdout=False, stderr=True)
extra_cp = snakemake.params.get("extra", "-v")
extra_iget = snakemake.params.get("extra_irods", "-vK")


# No node can access a cold storage
# these files must be copied. However,
# any file else where should be symlinked!
cold_storage = (
    "/mnt/isilon",
    "/mnt/archivage",
    "/mnt/install"
)


def cat_files(dest: str, *src: str, log: str = log) -> None:
    shell(f"cat {' '.join(src)} > {dest} {log}")


def bash_copy(src: str,
              dest: str,
              extra: str = extra_cp,
              log: str = log,
              cold: str = cold_storage) -> None:
    if src.startswith(cold) and ("--symbolic-link" not in extra):
        extra += " --symbolic-link"
    shell(f"cp {extra} {src} {dest} {log}")


def iRODS_copy(src: str,
               dest: str,
               extra: str = extra_iget,
               log: str = log) -> None:
    shell(f"iget {extra} {src} {dest} {log}")


def copy(src: str, dest: str) -> None:
    if src.startswith("/iRODS"):
        iRODS_copy(src, dest)
    else:
        bash_copy(src, dest)


def copy_then_concat(dest: str, *src: str) -> None:
    with TemporaryDirectory() as tmpdir:
        outfiles = []
        for path in src:
            copy(path, tmpdir)
            outfiles.append(op.join(tmpdir, op.basename(path)))
        cat_files(dest, *outfiles)


def copy_or_concat(src: str, dest: str) -> None:
    if "," in src:
        copy_then_concat(dest, *src.split(","))
    else:
        copy(src, dest)


output_directory = op.realpath(op.dirname(snakemake.output[0]))

sources = snakemake.params.get("input", snakemake.input)

if len(destinations := snakemake.output) == 1:
    # Then there is only one directory as a destination
    destination = op.realpath(str(destinations))
    shell(f"mkdir --parents --verbose {op.dirname(destination)}")
    if isinstance(sources, str):
        copy_or_concat(sources, destination)
    else:
        for source in sources:
            copy_or_concat(source, destination)
elif len(sources) == len(destinations):
    # Then there muse be as many paths as source
    for source, destination in zip(sources, destinations):
        copy_or_concat(source, destination)
else:
    print(len(sources), len(destinations), sources, destinations)
    raise ValueError("Number of input and output paths do not match.")
