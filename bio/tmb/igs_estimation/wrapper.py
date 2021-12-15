#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Compute the interrogated genomse seuqence size"""


import logging

from pybedtools import BedTool
from typing import Any, Dict, List, Union
from yaml import dump

logging.basicConfig(
    filename=snakemake.log[0],
    filemode="w",
    level=logging.DEBUG
)

def compute_feature_len(feature: List[Any]) -> List[Union[str, int]]:
    """
    Compute the length of a given feature
    """
    return (feature[0], int(feature[2]) - int(feature[1]))


def aggregate_per_chr(bedfile) -> Dict[str, Any]:
    logging.info("Computing per chromosome, genomic sequence size")
    lengths = iter(compute_feature_len(f) for f in bedfile)
    cchrom, clen = next(lengths, [None, None])
    current_chr = cchrom
    total_len = 0
    res_dict = {}
    while clen is not None:
        if current_chr == cchrom:
            total_len += clen
        else:
            logging.debug(f"New chromosome {current_chr} -> {cchrom}")
            res_dict[current_chr] = total_len
            total_len = clen
            current_chr = cchrom
        cchrom, clen = next(lengths, [None, None])
    res_dict[current_chr] = total_len
    res_dict["IGS"] = sum(res_dict.values())
    return res_dict



logging.info("Reading BED file")
bedfile = BedTool(snakemake.input["bed"]).merge()
igs = aggregate_per_chr(bedfile)
formatted_text = dump(igs, default_flow_style=False)

logging.info("Saving results")
with open(snakemake.output["yaml"], "w") as igs_out:
    igs_out.write(formatted_text)
