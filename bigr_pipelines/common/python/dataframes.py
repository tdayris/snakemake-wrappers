"""
This file contains usefull functions for pandas dataframe managment
"""

import itertools
import pandas
import logging

from typing import Optional, Union

def filter_df(df: pandas.DataFrame,
              colname: str,
              levels: list[str]) -> pandas.DataFrame:
    """
    Filter a dataframe given a column name and a list of levels to keep
    """
    return df[df[colname].isin(levels)]


def yield_comps(complete_design: pandas.DataFrame,
                aggregate: Optional[list[str]] = None,
                remove: Optional[list[str]] = None) \
                -> list[Union[pandas.DataFrame, list[str]]]:
    """
    Split a large design into simple ones and provide comparison information
    """
    if aggregate is not None:
        for cols in aggregate:
            complete_design["_".join(cols)] = (
                complete_design[cols].astype(str)
                                     .apply("_".join, axis=1)
                                     .str.strip()
                                     .str.replace(" ", "_")
            )

    if remove is not None:
        complete_design.drop(remove, axis=1, inplace=True)

    for col in complete_design.columns:
        levels = sorted(
            k for k, v in complete_design[col].value_counts().items() if v > 1
        )

        if len(levels) <= 1:
            # We need at least two levels to perfom a differential analysis
            continue


        for l1, l2 in itertools.combinations(levels, 2):
            # R and DESeq2 sort levels through alphanumerical order. We have
            # to provide user-understandable name ; so let us sort out the
            # levels and guess reference name.
            level, ref = sorted([l1, l2])

            # Building humand readable design name
            yield [col, level, ref]


def yield_samples(complete_design: pandas.DataFrame,
                  aggregate: Optional[list[str]] = None,
                  remove: Optional[list[str]] = None) \
                  -> list[Union[pandas.DataFrame, list[str]]]:
    """
    Split a large design into simple ones and provide comparison information
    """
    if aggregate is not None:
        for cols in aggregate:
            complete_design["_".join(cols)] = (
                complete_design[cols].astype(str)
                                     .apply("_".join, axis=1)
                                     .str.replace(" ", "_")
                                     .str.strip()
            )
    if remove is not None:
        complete_design.drop(remove, axis=1, inplace=True)

    for col in complete_design.columns:
        levels = sorted(
            k for k, v in complete_design[col].value_counts().items() if v > 1
        )

        if len(levels) <= 1:
            # We need at least two levels to perfom a differential analysis
            continue


        for l1, l2 in itertools.combinations(levels, 2):
            yield complete_design[complete_design[col].isin([l1, l2])].index.tolist()

def relation_condition_sample(complete_design: pandas.DataFrame,
                              factor: str) -> dict[str, str]:
    """
    From a design dataframe and a factor name, return the list of samples
    involved.
    """
    return dict(
        zip(complete_design.index.tolist(), complete_design[factor].tolist())
    )
