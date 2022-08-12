#!/usr/bin/env python3
# coding: utf-8

"""Write down material and methods for this pipeline"""

from typing import Optional

def is_default(params: Optional[str] = None) -> bool:
    if params is None:
        return True
    elif params.lower() == "none":
        return True
    
    return False

def parameters_explained(params: Optional[str] = None) -> str:
    if is_default(params) is True:
        return "default arguments"
    return "the following default parameters `{params}`"


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

Fastq file trimming was performeb by Fastp[1], using {fastp_parameters(config["fastp"]["adapters"])}
"""

with open(snakemake.output[0], "w") as mnm_stream:
    mnm_stream.write(text)