#!/usr/bin/env python
# coding: utf-8

import logging
import os
from symbol import comparison
import pandas
import seaborn
import matplotlib.pyplot

from pathlib import Path


def read_enrich(path: str) -> pandas.DataFrame:
    """Load an enrich result in memory"""
    comparison = os.sep.split(path)[-2]
    logging.info(f"Reading {path} in memory (cluster = {comparison}")
    tmp = pandas.read_csv(path, sep=" ", header=0, index_col=None)
    tmp["Comparison"] = [comparison for _ in tmp["p.adjust"]]
    return tmp


logging.basicConfig(filename=snakemake.log[0], filemode="w", level=logging.DEBUG)
seaborn.set_theme()

df = pandas.concat(read_enrich(path=p) for p in snakemake.input["tsv"])
df.sort_values(by="p.adjust", inplace=True)

seaborn.relplot(data=df, y="ID", x="p.adjust", hue="Comparison", size="Count")
matplotlib.pyplot.savefig(snakemake.output["relplot"], bbox_inches="tight")
matplotlib.pyplot.close()
