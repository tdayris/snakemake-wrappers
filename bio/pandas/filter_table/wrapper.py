#!/usr/bin/python3.8
# conding: utf-8

"""
Filter a TSV file
"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2020, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"


import logging
import operator
import pandas
import numpy
from typing import Union


def filter_dataframe(dataframe: pandas.DataFrame,
                     column: Union[str, int],
                     operator: str,
                     value: Union[int, float]):
    """
    Filter a dataframe according to the column, value and boolean
    operator
    """
    operator = ops[operator]
    return dataframe[operator(dataframe[column], value)]

logging.basicConfig(
    filename=snakemake.log[0],
    filemode="w",
    level=logging.DEBUG
)

ops = {
    ">": operator.gt,
    ">=": operator.ge,
    "==": operator.eq,
    "<": operator.lt,
    "<=": operator.le
}

if (sep := snakemake.params.get("separator", "\t")) != "\t":
    sep = snakemake.params["separator"]

data = pandas.read_csv(
    snakemake.input["table"],
    sep=sep,
    header=0
)

if (cols := snakemake.params.get("keep_column", None)) is not None:
    logging.debug(f"The following columns are kept: {cols}")
    data = data[cols]

if (line := snakemake.params.get("keep_line", None)) is not None:
    logging.debug(f"The following columns are kept: {line}")
    data = data[line]

if (filters := snakemake.params.get("filters", None)) is not None:
    logging.debug(f"The table will be filtered according to: {filters}")
    for filter in filters:
        data = filter_dataframe(data.copy, *filter)

logging.debug(f"Head of the final DataFrame:\n{data.head()}")

data.to_csv(
    snakemake.output["table"],
    sep=sep
)
