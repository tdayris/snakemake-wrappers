#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
This file contains useful functions for pandas data frame management
"""

import itertools
import pandas
import logging

from typing import Dict, List, Optional, Union


def filter_df(
    df: pandas.DataFrame, colname: str, levels: List[str]
) -> pandas.DataFrame:
    """
    Filter a dataframe given a column name and a list of levels to keep
    """
    return df[df[colname].isin(levels)]


def yield_comps(
    complete_design: pandas.DataFrame,
    aggregate: Optional[List[str]] = None,
    remove: Optional[List[str]] = None,
    contains: Optional[List[str]] = None,
) -> List[Union[pandas.DataFrame, List[str]]]:
    """
    Split a large design into simple ones and provide comparison information
    """
    if aggregate is not None:
        for cols in aggregate:
            complete_design["_".join(cols)] = (
                complete_design[cols]
                .astype(str)
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
            # We need at least two levels to perform a differential analysis
            continue

        for l1, l2 in itertools.combinations(levels, 2):
            # R and DESeq2 sort levels through alphanumerical order. We have
            # to provide user-understandable name ; so let us sort out the
            # levels and guess reference name.
            level, ref = sorted([l1, l2])
            for level, ref in [[l1, l2], [l2, l1]]:
                # Building human readable design name
                # Edit:
                # All previous attempt to guess reference name failed. This is now
                # me, not guessing anymore and doing all possible comparisons
                # two by two.
                try:
                    if any(
                        i in contains
                        for i in [
                            f"test_{level}",
                            f"reference_{ref}",
                            f"test_{l1}_vs_reference_{l2}",
                        ]
                    ):
                        yield [col, f"test_{level}", f"reference_{ref}"]
                except TypeError:
                    yield [col, f"test_{level}", f"reference_{ref}"]


def yield_samples(
    complete_design: pandas.DataFrame,
    aggregate: Optional[List[str]] = None,
    remove: Optional[List[str]] = None,
    contains: Optional[List[str]] = None,
) -> List[Union[pandas.DataFrame, List[str]]]:
    """
    Split a large design into simple ones and provide comparison information
    """
    if aggregate is not None:
        for cols in aggregate:
            complete_design["_".join(cols)] = (
                complete_design[cols]
                .astype(str)
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

        # print('aggregate:', aggregate)
        # print('current_column:', col)
        for l1, l2 in itertools.combinations(levels, 2):
            # Edit:
            # All previous attempt to guess reference name failed. This is now
            # me, not guessing anymore and doing all possible comparisons
            # two by two.
            try:
                if any(
                    i in contains
                    for i in [
                        f"test_{l1}",
                        f"reference_{l2}",
                        f"test_{l1}_vs_reference_{l2}",
                        f"test_{l2}_vs_reference_{l1}",
                    ]
                ):
                    yield complete_design[
                        complete_design[col].isin([l1, l2])
                    ].index.tolist()
                    # yield complete_design[complete_design[col].isin([l1, l2])].index.tolist()
            except TypeError:
                yield complete_design[
                    complete_design[col].isin([l1, l2])
                ].index.tolist()
                # yield complete_design[complete_design[col].isin([l1, l2])].index.tolist()


def relation_condition_sample(
    complete_design: pandas.DataFrame,
    factor: str,
    test: Optional[str] = None,
    ref: Optional[str] = None,
) -> Dict[str, str]:
    """
    From a design data frame and a factor name, return the list of samples
    involved.
    """
    if (test is None) and (ref is None):
        return dict(
            zip(complete_design.index.tolist(), complete_design[factor].tolist())
        )

    ref = ref[len("reference_") :] if ref.startswith("reference_") else ref
    test = test[len("test_") :] if test.startswith("test_") else test

    filtered_design = complete_design[complete_design[factor].isin([ref, test])]
    return dict(zip(filtered_design.index, filtered_design[factor]))
