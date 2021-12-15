#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Compute Tumor Mutation Burden given the definition of
Cancer Discov. 2020 Dec; 10(12): 1808–1825.
Published online 2020 Nov 2. doi: 10.1158/2159-8290.CD-20-0522

Tumor Mutational Burden (TMB) as a Predictive Biomarker in
Solid Tumors
"""

import logging
import os.path
import pandas
import yaml

logging.basicConfig(
    #filename=snakemake.log[0],
    filemode="w",
    level=logging.DEBUG
)


def build_sample_table(yamls: list[str]) -> pandas.DataFrame:
    sample_dict = {}
    for y in yamls:
        sample_name = os.path.basename(y)[:-len(".yaml")]
        with open(y, "r") as sample_yaml:
            data = yaml.safe_load(sample_yaml)
        sample_dict[sample_name] = data

    return pandas.DataFrame.from_dict(
        sample_dict,
        orient="index"
    )

logging.info("Loading datasets")
with open(snakemake.input["igs"], 'r') as igsyaml:
    igs = yaml.safe_load(igsyaml)
logging.debug(f"IGS: {igs}")
igs_mb = sum(igs.values()) / 1_000_000
logging.info(f"IGS size = {igs_mb}Mb")

samples = build_sample_table(snakemake.input["samples"])
logging.debug(f"Head of sample table:\n{samples.head()}")

logging.info("Computing TMB")
high_threshold = snakemake.params.get("high_threshold", 10)

for k, v in igs.items():
    try:
        samples[k] = (samples[k] / v) * 1_000_000
    except KeyError:
        samples[int(k)] = samples[int(k)] / igs_mb

samples["TMB"] = [
    f"{tmb} (TMB High)" if tmb > high_threshold
    else f"{tmb} (TMB Low)"
    for tmb in samples["IGS"]
]

samples["Mean"] = samples.drop(["IGS", "TMB"], axis=1).mean(axis=1)
samples["StDev"] = samples.drop(["IGS", "TMB"], axis=1).std(axis=1)

mean_chr = samples.drop(["TMB"], axis=1).mean(axis=0)
mean_chr.name = "Mean"

std_chr = samples.drop(["TMB"], axis=1).std(axis=0)
std_chr.name = "StdDev"
sample = samples.append([mean_chr, std_chr])

first_column = samples.pop("TMB")
samples.insert(0, "TMB", first_column)

logging.info("Saving results to disk")
samples.to_csv(snakemake.output["tsv"], sep="\t")
