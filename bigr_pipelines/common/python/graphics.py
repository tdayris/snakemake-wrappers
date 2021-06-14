"""
This file contains functions relative to graphs and image formats
"""

import math
import pandas

from typing import Optional

def optimal_size(multiplier: int = 1) -> list[int]:
    """
    Return optimal width/height of a graph
    """
    return [1024 * multiplier, 768 * multiplier]


def image_size_from_sample_number(sample_nb: Optional[int] = None) -> list[int]:
    """
    Return usual image size given a sample number
    """
    return optimal_size(multiplier=min(math.ceil(sample_nb / 4), 1))


def image_size_from_design(design: pandas.DataFrame,
                           factor: str,
                           test: str,
                           ref: str) -> list[int]:
    """
    Return usual image size given a design and a comparison level
    """
    return image_size_from_sample_number(
        sample_nb=len(design[design[factor].isin([test, ref])][factor].tolist())
    )
