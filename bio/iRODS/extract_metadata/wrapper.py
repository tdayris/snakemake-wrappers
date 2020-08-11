#!/usr/bin/python3.7
# -*- coding: utf-8 -*-

"""Snakemake wrapper for iRODS metadata extraction"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2020, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

import pandas
import yaml

from typing import Any, Dict, List, Union

# Load complete metadata file
with open(snakemake.input[0], "r") as yamlin:
    metadatas = yaml.load(yamlin, Loader=yaml.SafeLoader)


# Extract user defined metadatas
def getitem(item: Union[str, List[str]],
            indict: Dict[str, Any]) -> Any:
    """
    Return a list of keys from a dictionnary, providing None as default
    missing value
    """
    if isinstance(item, str):
        try:
            return {item: indict[item]['value']}
        except KeyError:
            return None

    results = {}
    for i in item:
        try:
            results[i] = indict[i]["value"]
        except KeyError:
            results[i] = None

    return results


content = [
    {"Path": k, **getitem(snakemake.params["attributes"], v)}
    for k, v in metadatas.items()
]

# Formatting output as TSV, then saving it.
pandas.DataFrame(content).to_csv(
    snakemake.output[0],
    sep="\t",
    header=True,
    index=True
)
