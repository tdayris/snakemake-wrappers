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
import re
from typing import Any, Optional, Union


def filter_dataframe(dataframe: pandas.DataFrame,
                     column: Union[str, int],
                     local_operator: str,
                     value: Union[int, float]) -> pandas.DataFrame:
    """
    Filter a dataframe according to the column, value and boolean
    operator
    """
    logging.info(
        "Filtering %s with operator %s", str(column), str(local_operator)
    )
    try:
        local_operator = ops[local_operator]
        return dataframe[local_operator(dataframe[column], value)]
    except TypeError:
        logging.error(dataframe[column].head())
        raise



def filter_full_lines(dataframe: pandas.DataFrame,
                      local_operator: str,
                      value: Union[int, float]) -> pandas.DataFrame:
    """
    Apply filters on whole lines
    """
    local_operator = ops[local_operator]
    return data.loc[~local_operator(data, value).all(axis=1)]


def add_column(dataframe: pandas.DataFrame,
               new_column: Union[str, int],
               local_operator: str,
               input_col: Union[str, int],
               value: Optional[Union[str, int]] = None) -> pandas.DataFrame:
    """
    Create a new column based on information from the first one, and the
    provided operator and value.
    """
    logging.debug("Operator was: %s", local_operator)
    logging.debug("On: %s", str(input_col))
    if local_operator in ["prefix", "suffix"]:
        dataframe[input_col] = dataframe[input_col].astype(str)
    local_operator = ops[local_operator]
    if value is None:
        dataframe[new_column] = local_operator(dataframe[input_col])
    else:
        logging.debug("Value is: %s", str(value))
        dataframe[new_column] = local_operator(dataframe[input_col], value)
    return dataframe


def combine_columns(dataframe: pandas.DataFrame,
                    new_column: Union[str, int],
                    local_operator: str,
                    left_input_col: Union[str, int],
                    right_input_col: Union[str, int]) -> pandas.DataFrame:
    """
    Create a new column based on information from the first one, and the
    provided operator, left and right columns.
    """
    local_operator = ops[local_operator]
    dataframe[new_column] = local_operator(
        dataframe[left_input_col],
        dataframe[right_input_col]
    )
    return dataframe


def prepare_int_conversion(value: Any) -> Any:
    if isinstance(value, (int, float)) or value.isnumeric():
        return value
    if value.upper() in ["", ".", "NA", "NAN", "#"]:
        logging.warning(f"Could not convert {value} to int/float")
        return numpy.nan
    try:
        value = value.strip(",").strip("|")
        if "|" in value:
            logging.warning(f"Multiple values in {value}, the last one has been kept")
            return float([i for i in value.split("|") if i != ""][-1])
        logging.warning(f"Multiple values in {value}, the last one has been kept")
        return float([i for i in value.split(",") if i != ""][-1])
    except ValueError:
        logging.error(f"Failed converting string to float {value}")
        return numpy.nan



logging.basicConfig(
    filename=snakemake.log[0],
    filemode="w",
    level=logging.DEBUG
)


ops = {
    ">": pandas.Series.gt,
    ">=": pandas.Series.ge,
    "==": pandas.Series.eq,
    "<": pandas.Series.lt,
    "<=": pandas.Series.le,
    "!=": pandas.Series.ne,
    "+": pandas.Series.add,
    "-": pandas.Series.sub,
    "/": pandas.Series.div,
    "*": pandas.Series.mul,
    "//": pandas.Series.floordiv,
    "abs": pandas.Series.abs,
    "prefix": pandas.Series.add_prefix,
    "suffix": pandas.Series.add_suffix,
    "as": pandas.Series.astype,
    "=": pandas.Series.copy
}


if (sep := snakemake.params.get("separator", "\t")) != "\t":
    sep = snakemake.params["separator"]


data = pandas.read_csv(
    snakemake.input["table"],
    sep=sep,
    header=0
)


if (new_cols := snakemake.params.get("new_cols", None)) is not None:
    logging.debug("Creating the following new columns: %s", str(new_cols))
    for new_col in new_cols:
        data = add_column(data, *new_col)

if (combine_cols := snakemake.params.get("combine_cols", None)) is not None:
    logging.debug("Creating a the following new columns, as a combination of existing columns: %s", str(combine_cols))
    for combine_col in combine_cols:
        data = combine_column(data, *combine_col)


if (prefixes := snakemake.params.get("prefixes", None)) is not None:
    logging.debug("Adding the following prefixes %s", str(prefixes))
    for col, prefix in prefixes:
        data[col] = [f"{prefix}{val}" for val in data[col].astype(str)]


if (cols := snakemake.params.get("keep_column", None)) is not None:
    logging.debug(f"The following columns are kept: {cols}")
    data = data[cols]


if (not_cols := snakemake.params.get("drop_column", None)) is not None:
    logging.debug(f"The following columns are dropped out: {not_cols}")
    data = data[list(set(data.columns.tolist()) - set(not_cols))]


if (convert_cols_type := snakemake.params.get("convert_cols_type", None)) is not None:
    logging.debug(f"The following type concersion are made: {convert_cols_type}")
    for column_name, new_type in convert_cols_type.items():
        if new_type in ["int", "float"]:
            logging.info(f"Converting {column_name}")
            data[column_name] = [
                prepare_int_conversion(i) for i in data[column_name]
            ]
        data[column_name] = data[column_name].astype(new_type)


if (line := snakemake.params.get("keep_line", None)) is not None:
    logging.debug(f"The following lines are kept: {line}")
    data = data[data.index.isin(line)]


if (not_line := snakemake.params.get("drop_line", None)) is not None:
    logging.debug(f"The following lines are droped out {not_line}")
    data = data[not data.index.isin(not_line)]


if (filters := snakemake.params.get("filters", None)) is not None:
    for filter in filters:
        logging.debug(f"The table will be filtered according to: {filter}")
        data = filter_dataframe(data, *filter)


if (filters := snakemake.params.get("full_line_filters", None)) is not None:
    for filter in filters:
        logging.debug(f"Applying the following filter on whole lines: {filter}")
        data = filter_full_lines(data, *filter)


if (not_contains := snakemake.params.get("not_contains", None)) is not None:
    for column, value in not_contains:
        logging.debug(f"Filtering out line in which {column} contains: {value}")
        data = data[~data[column].str.contains(value, flags=re.IGNORECASE, regex=True)]


if (contains := snakemake.params.get("contains", None)) is not None:
    for column, value in contains:
        logging.debug(f"Filtering in line in which {column} contains: {value}")
        data = data[data[column].str.contains(value, flags=re.IGNORECASE, regex=True)]


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
    logging.debug(f"Setting {new_index} as first column")
    index_values = data.pop(new_index)
    data.insert(new_index, 0, index_values, allow_duplicates = False)


logging.debug(f"Head of the final DataFrame:\n{data.head()}")


if (new_index := snakemake.params.get("set_index", None)) is not None:
    previous_index = None
    if snakemake.params.get("drop_duplicated_index", True) is True:
        data.drop_duplicates(subset=new_index, keep="first", inplace=True)

    if snakemake.params.get("override_previous_index", False) is True:
        previous_index = "index" if data.index.name is None else str(data.index.name)

    data.reset_index(inplace=True)
    data.set_index(new_index, inplace=True)

    try:
        del data[previous_index]
    except KeyError:
        pass

    if (not_cols := snakemake.params.get("drop_column", None)) is not None:
        logging.debug(f"The following columns are dropped out: {not_cols}")
        data = data[list(set(data.columns.tolist()) - set(not_cols))]


data.to_csv(
    snakemake.output["table"],
    sep=sep,
    index=snakemake.params.get("keep_index", False)
)

if "xlsx" in snakemake.output.keys():
    data.to_excel(
        snakemake.output["xlsx"],
        header=True,
        index=True,
        na_rep="."
    )
