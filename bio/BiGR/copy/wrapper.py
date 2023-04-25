#!/usr/bin/python3.7
# -*- coding: utf-8 -*-

"""Snakemake wrapper for bash copy within the IGR's Flamingo cluster"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2019, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

import os
import os.path as op
import logging

from snakemake.shell import shell
from tempfile import TemporaryDirectory

# Logging behaviour
try:
    logging.basicConfig(
        filename=snakemake.log[0],
        filemode="w",
        level=logging.DEBUG
    )
except Exception:
    logging.basicConfig(level=logging.DEBUG)
    logging.warning(
        "No logging file was provided in Snakemake rule. "
        "Logging will be written in standard output."
    )


# Prepare logging
log = snakemake.log_fmt_shell(stdout=False, stderr=True, append=True)
extra_cp = snakemake.params.get("extra", "--force --recursive --verbose")
extra_iget = snakemake.params.get("extra_irods", "-vKrf")
extra_ln = snakemake.params.get("extra_ln", "--relative --verbose --symbolic --force")
extra_rsync = snakemake.params.get(
    "extra_rsync", 
    "--verbose --recursive --checksum --human-readable --partial --progress"
)
use_cp_over_rsync = snakemake.params.get("use_cp_over_rsync", False)
master_node = snakemake.params.get("master_node", "flamingo-lg-01")
parallel = snakemake.params.get("internal_parallel", False)

if ("-N" not in extra_iget) and snakemake.threads > 1:
    max_threads = snakemake.threads
    if os.uname()[1] == master_node and max_threads > 4:
        logging.warning(
            "More than 4 threads per copy on master node is not allowed. "
            "Max threads set to 4. This does not change your original "
            "threading reservation. Run this script on a computing node, "
            "or lower the number of threads."
        )
        max_threads = 4
    if parallel:
        logging.info("Internal parallelization")
        max_threads = 1
    extra_iget += f" -N {max_threads} "

# No node can access a cold storage
# these files must be copied. However,
# any file else where should be symlinked!
cold_storage = (
    "/mnt/isilon",
    "/mnt/archivage",
    "/mnt/install"
)
if "cold_storage" in snakemake.params.keys():
    cold_storage = snakemake.params["cold_storage"]


def cat_files(dest: str, *src: str, log: str = log) -> None:
    command = f"cat {' '.join(src)} > {dest} {log}"
    logging.info(command)
    shell(command)


def bash_copy(src: str,
              dest: str,
              extra_cp: str = extra_cp,
              extra_rsync: str = extra_rsync,
              extra_ln: str = extra_ln,
              log: str = log,
              cold: str = cold_storage,
              use_cp: bool = use_cp_over_rsync) -> None:
    if not src.startswith(cold):
        command = f"ln {extra_ln} {src} {dest} {log}"
        logging.info(command)
        shell(command)
    else:
        command ='cp' if use_cp else 'rsync'
        extra = extra_cp if use_cp else extra_rsync
        command = f"{command} {extra} {src} {dest} {log}"
        logging.info(command)
        shell(command)


def iRODS_copy(src: str,
               dest: str,
               extra: str = extra_iget,
               log: str = log) -> None:
    command = f"iget {extra} {src} {dest} {log}"
    logging.info(command)
    shell(command)



def copy(src: str, dest: str) -> None:
    if src.startswith("/odin/kdi/dataset/"):
        iRODS_copy(src, dest)
    else:
        bash_copy(src, dest)


def copy_then_concat(dest: str, *src: str) -> None:
    with TemporaryDirectory() as tmpdir:
        outfiles = []
        for path in src:
            tmp_dest = f"{tmpdir}/{os.path.basename(path)}.{os.urandom(8)}"
            copy(path, tmp_dest)
            outfiles.append(tmp_dest)
        cat_files(dest, *outfiles)


def copy_or_concat(src: str, dest: str) -> None:
    if "," in src:
        copy_then_concat(dest, *src.split(","))
    else:
        copy(src, dest)


output_directory = op.realpath(op.dirname(snakemake.output[0]))

sources = snakemake.params.get("input", snakemake.input)
destinations = snakemake.output

if len(destinations) == 1:
    # Then there is only one directory as a destination
    destination = op.realpath(str(destinations))
    command = f"mkdir --parents --verbose {op.dirname(destination)} {log}"
    logging.info(command)
    shell(command)
    if isinstance(sources, str):
        copy_or_concat(sources, destination)
    else:
        for source in sources:
            copy_or_concat(source, destination)
elif len(sources) == len(destinations):
    # Then there must be as many paths as source
    for source, destination in zip(sources, destinations):
        command = f"mkdir --parents --verbose {op.dirname(destination)} {log}"
        logging.info(command)
        shell(command)
        copy_or_concat(source, destination)
else:
    logging.error(f"nb of input: {len(sources)}, nb of output: {len(destinations)}")
    logging.debug(sources)
    logging.debug(destinations)
    raise ValueError("Number of input and output paths do not match.")
