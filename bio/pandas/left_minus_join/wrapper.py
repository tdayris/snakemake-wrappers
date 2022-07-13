#!/usr/bin/python3.8
# conding: utf-8

"""
Perform a left minus join on two dataframes with python
"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2020, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"


import pandas
import logging

from typing import List, Optional, Union

logging.basicConfig(
    filename=snakemake.log[0],
    filemode="w",
    level=logging.DEBUG
)

def read_table(path: str,
               separator: str = "\t",
               index_cols: Optional[Union[str, List[str], int]] = 0,
               header_col: Optional[Union[str, int]] = None) -> pandas.DataFrame:
    """
    Read a provided table and return an indexed dataframe
    """
    data = pandas.read_csv(
        path,
        sep=separator,
        header=header_col
    )
    data.set_index(
        index_cols,
        inplace=True
    )
    logging.debug("Head of dataframe:")
    logging.debug(data.head())
    return data

# Loading optional parameters
sep = snakemake.params.get("sep", "\t")
index = snakemake.params.get("index", 0)
header = snakemake.params.get("header", None)
logging.debug(
    "Optional parameters: separator = {}, index = {}, header = {}".format(
        sep, index, header
    )
)

# Loadin input datasets
left_table = read_table(
    path=snakemake.input["left"],
    separator=sep,
    header_col=header,
    index_cols=index
)

right_table = read_table(
    path=snakemake.input["right"],
    separator=sep,
    header_col=header,
    index_cols=index
)

# Perfoming left minus join
left_minus_join = left_table[
    left_table.index.isin(right_table.index) == False
]

# Saving results
left_minus_join.to_csv(
    snakemake.output[0],
    sep=sep,
    header=True,
    index=True
)
logging.debug("Head of join result:")
logging.debug(left_minus_join.head())
