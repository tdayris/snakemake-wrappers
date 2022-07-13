#!/usr/bin/env python3
# coding: utf-8

"""Merge multiple TSV with pandas"""

import logging
import pandas

from typing import Any

try:
    # Case user provided logging file
    logging.basicConfig(
        filename=snakemake.log[0],
        filemode="w",
        level=logging.DEBUG
    )
except IndexError:
    # Case user did not provide any logging file
    logging.basicConfig(filemode="w", level=logging.DEBUG)

    
def read_single_tsv(path: str, extra: dict[str, Any]) -> pandas.DataFrame:
    """Return a pandas.DataFrame from a path to a TSV file"""
    logging.debug("Loading %s", path)
    return pandas.read_csv(
        path,
        **extra
    )


def merge_tsv(paths: list[str], 
              extra_read: dict[str, Any], 
              extra_merge: dict[str, Any]
             ) -> pandas.DataFrame:
    """Read multiple TSV files and merge them according to the extra params"""
    result = None
    for path in paths:
        tmp = read_single_tsv(path, extra_read)
        try:
            result = pandas.merge(
                left=result,
                right=tmp,
                **extra_merge
            )
        except AttributeError:
            # The first and only the first merging attempt will fail,
            # because result is None. Else, there should not be any issue.
            retult = tmp

    logging.debug("Tables merged: %s", str(result.head()))
    return result


if len(snakemake.input["tsv"]) <= 1:
    raise ValueError("Cannot merge 1 or less file.")

df = merge_tsv(
    paths=snakemake.input["tsv"],
    extra_read=snakemake.params.get("extra_read", {}),
    extra_merge=snakemake.params.get("extra_merge", {})
)

logging.debug("Saving results to disk")
df.to_csv(
    snakemake.output[0],
    sep="\t",
    **snakemake.params.get("extra_write", {})
)