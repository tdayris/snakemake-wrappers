#!/usr/bin/python3.9
# -*- coding: utf-8 -*-

"""This wrapper aggregates multiple .msi files from msi-sensor"""

import pandas

def read_msi(path: str, sample: str) -> pandas.DataFrame:
    tmp = pandas.read_csv(path, index_col=None, header=0, sep="\t")
    tmp["Sample_id"] = [sample]
    tmp.set_index("Sample_id", inplace=True)
    return tmp


def stability(percent: float, stability: float) -> list[str]:
    if percent >= stability:
        return "MSI"
    return "MSS"


msi = None
for msi_path, sample_name in zip(snakemake.input, snakemake.params.sample_list):
    tmp = read_msi(msi_path, sample_name)
    try:
        msi = pandas.concat([msi, tmp], axis=0, verify_integrity=True)
    except ValueError:
        msi = tmp

msi["Stability"] = [
    stability(percent, snakemake.params.get("stability", 20)) 
    for percent in msi["%"]
]

msi.to_csv(
    snakemake.output[0],
    sep="\t",
    header=True,
    index=True,
)