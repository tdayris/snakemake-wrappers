#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Snakemake wrapper for subsampling"""

from functools import singledispatch
from random import uniform
from typing import Dict, List

from gzip import open as gzopen

@singledispatch
def fraction(expected: float, total: float) -> float:
    """From an expected number of reads, return input file fraction"""
    return expected / total

@fraction.register(list)
def fraction(expected: List[float], total: float) -> float:
    """From a mean and std, return input file fraction"""
    vmin = expected[0] - expected[1]
    vmax = expected[0] + expected[1]
    return random.uniform(vmin, vmax) / total


@fraction.register(dict)
def fraction(expected: Dict[str, float], total: float) -> float:
    """From a mean and std, return input file fraction"""
    vmin = expected["mean"] - expected["std"]
    vmax = expected["mean"] + expected["std"]
    return random.uniform(vmin, vmax) / total


def nb_reads(fq_path: str) -> float:
    with gzopen(fq_path) as fq_stream:
        return sum(1 for _ in fq_stream) / 4


total = nb_reads(snakemake.input["fastq"])
subsample = fraction(
    snakemake.params.get("read_nb", {"mean": 6.25576e+07, "std": 1.05253e+07})
    total
)

try:
    seed = snakemake.wildcards["seed"]
except KeyError:
    seed = snakemake.params.get("seed", random.randint(0, 999))

shell(
    "sambamba view --nthreads {snakemake.threads} --subsample {subsample} "
    "--subsample-seed {seed} "
)
