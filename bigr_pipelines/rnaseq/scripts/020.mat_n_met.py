#!/usr/bin/env python3
# coding: utf-8

"""Write down material and methods for this pipeline"""

from typing import Optional

def is_default(params: Optional[str] = None) -> bool:
    if params is None:
        return True
    elif params.lower() in ["none", "null", "None", "NULL"]:
        return True
    
    return False

def parameters_explained(params: Optional[str] = None) -> str:
    if is_default(params) is True:
        return "default arguments"
    return "the following default parameters: `{params}`"


def fastp_parameters(adapters: Optional[str] = None, extra: Optional[str] = None) -> str:
    result = ""
    if not is_default(adapters):
        result = "custom adapter list, and using "
    result += parameters_explained(extra)

    return result

config = snakemake.params["all_parameters"]
text = f"""
# Material and Methods

## Quality controls

Fastq file trimming and raw quality control was performed by Fastp[1] (version: 0.20), using {fastp_parameters(config["fastp"]["adapters"])}. 
FastqScreen[2] (version: 0.5.2) was ran through trimmed fastq files to assess their specificity to exptected genomes (here: {config["reference"]["genome_name"]}).
Through the whole pipeline, quality control reports were aggregated in a single HTML file using MultiQC[3] (version: 1.10.1) using default parameters, and an in-house header. Tools and graphs that are not available natively on MultiQC were included through configuration file and custom graphs method, using in-house scripts.

## DGE

## Fusions

## MSI

## Deconvolution


## Citations

1. Shifu Chen, Yanqing Zhou, Yaru Chen, Jia Gu; fastp: an ultra-fast all-in-one FASTQ preprocessor, Bioinformatics, Volume 34, Issue 17, 1 September 2018, Pages i884–i890, https://doi.org/10.1093/bioinformatics/bty560
2. FastqScreen: Wingett SW, Andrews S. FastQ Screen: A tool for multi-genome mapping and quality control. F1000Res. 2018 Aug 24;7:1338. doi: 10.12688/f1000research.15931.2. PMID: 30254741; PMCID: PMC6124377.
3. Philip Ewels, Måns Magnusson, Sverker Lundin and Max Käller MultiQC: Summarize analysis results for multiple tools and samples in a single report. Bioinformatics (2016) doi: 10.1093/bioinformatics/btw354 PMID: 27312411 

"""

with open(snakemake.output[0], "w") as mnm_stream:
    mnm_stream.write(text)