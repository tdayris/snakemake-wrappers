#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Filter a VCF file"""

import logging
import gzip
import yaml

logging.basicConfig(
    filename=snakemake.log[0],
    filemode="w",
    level=logging.DEBUG
)

filter_in = set(snakemake.params.get("filter_in", []))
apply_filters_in = len(filter_in) > 0
filter_out = set(snakemake.params.get("filter_out", []))
apply_filters_out = len(filter_out) > 0

columns = None
tumor_sample = None
normal_sample = None

min_coverage = snakemake.params.get("min_coverage", 10)
format_corevage = snakemake.params.get("allele_depth", "AD")

result = {}

with gzip.open(snakemake.input["vcf"], "r") as vcfin:
    for line in vcfin:
        line = line.decode("utf-8")

        if line.startswith("##"):
            continue

        if line.startswith("#"):
            columns = line[1:-1].split("\t")
            tumor_sample = columns[-1]
            logging.info(f"Tumor sample: {tumor_sample}")
            normal_sample = columns[-2]
            logging.info(f"Normal sample: {normal_sample}")
            continue

        chomp = dict(zip(columns, line[:-1].split("\t")))
        if chomp["CHROM"] not in result.keys():
            result[chomp["CHROM"]] = 0

        filters = set(chomp["FILTER"].split(","))

        if apply_filters_out:
            if len(filters & filter_out) > 0:
                logging.debug(f"{chomp['CHROM']}:{chomp['POS']}-{chomp['REF']}>{chomp['ALT']} filtered out on {chomp['FILTER']}")
                continue

        if apply_filters_in:
            if len(fitlers & filter_in) == 0:
                logging.debug(f"{chomp['CHROM']}:{chomp['POS']}-{chomp['REF']}>{chomp['ALT']} not filtered in on {chomp['FILTER']}")
                continue

        coverage = chomp["FORMAT"].split(":").index(format_corevage)
        tumor_coverage = chomp[tumor_sample].split(":")[coverage].split(",")[1]
        if tumor_coverage == ".":
            continue
        elif int(tumor_coverage) < min_coverage:
            logging.debug(f"{chomp['CHROM']}:{chomp['POS']}-{chomp['REF']}>{chomp['ALT']} not enough covered: {tumor_coverage}")
            continue

        normal_coverage = chomp[normal_sample].split(":")[coverage].split(",")[1]
        if normal_coverage != ".":
            if int(normal_coverage) >= min_coverage:
                logging.debug(f"{chomp['CHROM']}:{chomp['POS']}-{chomp['REF']}>{chomp['ALT']} was germline: {normal_coverage}")
                continue

        result[chomp["CHROM"]] += 1

result["IGS"] = sum(result.values())
formatted_text = yaml.dump(result, default_flow_style=False)

with open(snakemake.output["yaml"], "w") as igs_out:
    igs_out.write(formatted_text)
