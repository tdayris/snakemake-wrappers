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
                     local_operator: str,
                     value: Union[int, float]):
    """
    Filter a dataframe according to the column, value and boolean
    operator
    """
    local_operator = ops[local_operator]
    return dataframe[local_operator(dataframe[column], value)]


def filter_full_lines(dataframe: pandas.DataFrame,
                      local_operator: str,
                      value: Union[int, float]):
    """
    Apply filters on whole lines
    """
    local_operator = ops[local_operator]
    return data.loc[~local_operator(data, value).all(axis=1)]


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
    "<=": operator.le,
    "!=": operator.ne
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


if (not_cols := snakemake.params.get("drop_column", None)) is not None:
    logging.debug(f"The following columns are dropped out: {not_cols}")
    data = data[list(set(data.columns.tolist()) - set(not_cols))]


if (line := snakemake.params.get("keep_line", None)) is not None:
    logging.debug(f"The following lines are kept: {line}")
    data = data[data.index.isin(line)]


if (not_line := snakemake.params.get("drop_line", None)) is not None:
    logging.debug(f"The following lines are droped out {not_line}")
    data = data[not data.index.isin(not_line)]


if (filters := snakemake.params.get("filters", None)) is not None:
    logging.debug(f"The table will be filtered according to: {filters}")
    for filter in filters:
        data = filter_dataframe(data.copy, *filter)


if (filters := snakemake.params.get("full_line_filters", None)) is not None:
    logging.debug(f"Applying the following filter on whole lines: {filters}")
    for filter in filters:
        data = filter_full_lines(data.dopy, *filter)


if snakemake.params.get("dropna", False) is True:
    logging.debug("Dropping lines with NAs")
    data.dropna(axis=0, how="any", inplace=True)


if snakemake.params.get("drop_null_lines", False) is True:
    logging.debug("Dropping lines full of zeros")
    data = data.loc[~(data==0).all(axis=1)]


if snakemake.params.get("drop_na_lines", False) is True:
    logging.debug("Dropping lines full of NAs")
    data.dropna(axis=0, how="all", inplace=True)


if snakemake.params.get("drop_duplicated_lines", False) is True:
    logging.debug("Dropping duplicated lines")
    data.drop_duplicates(inplace=True)
    

if (dedup_cols := snakemake.params.get("drop_duplicates_on", None)) is not None:
    logging.debug(
        "Dropping duplicates on the following column(s): %s", str(dedup_cols)
    )
    data.drop_duplicates(subset=dedup_cols, keep="first", inplace=True)


if (new_index := snakemake.params.get("new_index_col", None)) is not None:
    logging.debug("Setting {nex_index} as first column")
    index_values = data.pop(new_index)
    data.insert(new_index, 0, index_values, allow_duplicates = False)


logging.debug(f"Head of the final DataFrame:\n{data.head()}")


if ("set_index" in snakemake.params.keys()):
    (data.reset_index(inplace=True)
         .set_index(snakemake.params["set_index"], inplace=True))

data.to_csv(
    snakemake.output["table"],
    sep=sep,
    index=snakemake.params.get("keep_index", False)
)
