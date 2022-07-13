#!/usr/bin/python3.9
# -*- coding: utf-8 -*-

"""
Compute md5hashes of provided files using a pool of actors.

Example:
python3 ./md5compute.py -t 5 *.fq.gz
"""

import hashlib
import io
import logging
import os
import pprint
import pykka
import sys


class Md5Computer(pykka.ThreadingActor):
    """
    The actor designed to compute md5 hashes of a given file
    """
    def on_start(self) -> None:
        """
        Logging method automatically called on start
        """
        logging.debug("%s: Starting a Md5Computer actor", self)

    def on_stop(self) -> None:
        """
        Logging method automatically called on stop
        """
        logging.debug("%s: Stopping a Md5Computer actor", self)

    def md5(self, path: str, buffer_length=io.DEFAULT_BUFFER_SIZE) -> list[str]:
        """
        Read files through a chunk to avoid too large file being stored into
        memory, and compute a md5sum of this file.
        """
        logging.debug("%s: Hashing %s", self, path)
        md5 = hashlib.md5()
        with io.open(path, mode="rb") as fd:
            for chunk in iter(lambda: fd.read(buffer_length), b''):
                md5.update(chunk)
        return md5.hexdigest()


def run(pool_size: int, *paths: str) -> list[list[str]]:
    """
    Creates the pool of actors, distribute the work, and gather results.
    """
    computers = [Md5Computer.start().proxy() for _ in range(pool_size)]
    hosts = [
        computers[idx % len(computers)].md5(path)
        for idx, path in enumerate(paths)
    ]
    hashes = zip(paths, pykka.get_all(hosts))
    pprint.pprint(list(hashes))
    pykka.ActorRegistry.stop_all()

    return hashes


def save_hashes(output_path: str, *hashes: list[str]) -> None:
    """
    Save results to a text file
    """
    with open(output_path, "w") as outfile:
        for file, hex in hashes:
            outfile.write(f"{hex}\t{os.path.basename(file)}\n")


# Python loggings
verbose = (
    logging.DEBUG
    if snakemake.params.get("verbose", False)
    else logging.INFO
)

level_logging = (logging.DEBUG if verbose else logging.INFO)
logging.basicConfig(
    filename=snakemake.log[0],
    filemode="w",
    level=level_logging
)
logging.getLogger('pykka').setLevel(logging.WARNING)

try:
    hexdigest_list = run(snakemake.threads, *snakemake.input)
    save_hashes(snakemake.output[0], *hexdigest_list)
    if any(hashval > 1 for hashval in Counter(hexdigest_list).values()):
        error_text = "Some input are duplicated"
        logging.error(error_text)
        raise ValueError(error_text)
except Exception as e:
    logging.error(e)
    raise
