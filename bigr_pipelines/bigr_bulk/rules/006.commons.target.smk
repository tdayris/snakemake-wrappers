"""
This snakefile contains python functions for:
* Listing possible target files
* Returning target files for Snakemake
"""

import functools
import pandas

from typing import Any, Dict, List

def contains_protocol(design: pandas.DataFrame, protocol: List[str]) -> bool:
    """
    Return True if design contains data belonging to the given protocol
    """
    return any(i for i in design["Prococol"].str.lower.isin(protocol))

contains_rnaseq = functools.partial(protocol=["rnaseq", "rna-seq"])
contains_wes = functools.partial(protocol=["wes", "exome", "exome-seq"])

def get_steps(config: Dict[str, Any] = config) -> bool:
    """
    Return list of steps from config and define default ones
    in case of missing keys
    """
    results = []
    steps = config.get("steps", {})
    if steps.get("Index", False):
        results.append("Index")
    if steps.get("Quantification", False):
        results.append("Quantification")
    
    return results

def get_targets() -> Dict[str, str]:
    """
    Return a dictionary of expected output files at the end
    of the BiGR bulk pipeline.
    """
    targets = {
        "trimming_qc": "data_output/MultiQC/Trimming.html"
    }

    if contains_rnaseq():
        pass

    if contains_wes():
        pass


    return targets