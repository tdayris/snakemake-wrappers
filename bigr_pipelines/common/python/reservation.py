#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import functools


# Memory and time reservation
def get_resources_per_gb(
    wildcards, input, attempt, multiplier: int = 1, base: int = 0
) -> int:
    """
    Return the amount of resources needed per GB of input.

    Parameters:
    * wildcards: Snakemake signature requires this parameter.
    * input: The input file list
    * attempt: The # of times the calling rule has been restarted
    * multiplier: An arbitrary multiplier

    Return:
    (int) The amount of resources needed (mb, minutes, etc)
    """
    return max(
        # Case there is 1gb or more in input
        # ((input.size // 10_000_000_000) * attempt * multiplier) + base,
        1,
        # Case there is less than 1gb in input
        (multiplier * attempt) + base,
    )


# Explicit time reservations
get_15min_per_attempt = functools.partial(get_resources_per_gb, multiplier=15)
get_20min_per_attempt = functools.partial(get_resources_per_gb, multiplier=20)
get_35min_per_attempt = functools.partial(get_resources_per_gb, multiplier=35)
get_45min_per_attempt = functools.partial(get_resources_per_gb, multiplier=45)
get_75min_per_attempt = functools.partial(get_resources_per_gb, multiplier=75)
get_1h_per_attempt = functools.partial(get_resources_per_gb, multiplier=60)
get_2h_per_attempt = functools.partial(get_resources_per_gb, multiplier=60 * 2)
get_3h_per_attempt = functools.partial(get_resources_per_gb, multiplier=60 * 3)
get_4h_per_attempt = functools.partial(get_resources_per_gb, multiplier=60 * 4)
get_5h_per_attempt = functools.partial(get_resources_per_gb, multiplier=60 * 5)
get_6h_per_attempt = functools.partial(get_resources_per_gb, multiplier=60 * 6)
get_8h_per_attempt = functools.partial(get_resources_per_gb, multiplier=60 * 8)
get_10h_per_attempt = functools.partial(get_resources_per_gb, multiplier=60 * 10)
get_16h_per_attempt = functools.partial(get_resources_per_gb, multiplier=60 * 16)
get_90min_per_attempt = functools.partial(get_resources_per_gb, multiplier=90)

# Explicit memory reservation
get_1gb_per_attempt = functools.partial(get_resources_per_gb, multiplier=1024)
get_2gb_per_attempt = functools.partial(get_resources_per_gb, multiplier=1024 * 2)
get_4gb_per_attempt = functools.partial(get_resources_per_gb, multiplier=1024 * 4)
get_5gb_per_attempt = functools.partial(get_resources_per_gb, multiplier=1024 * 5)
get_6gb_per_attempt = functools.partial(get_resources_per_gb, multiplier=1024 * 6)
get_8gb_per_attempt = functools.partial(get_resources_per_gb, multiplier=1024 * 8)
get_10gb_per_attempt = functools.partial(get_resources_per_gb, multiplier=1024 * 10)
get_15gb_per_attempt = functools.partial(get_resources_per_gb, multiplier=1024 * 15)
get_20gb_per_attempt = functools.partial(get_resources_per_gb, multiplier=1024 * 20)
get_75gb_and_2gb_per_attempt = functools.partial(
    get_resources_per_gb, multiplier=1024 * 2, base=1024 * 75
)
get_75gb_and_5gb_per_attempt = functools.partial(
    get_resources_per_gb, multiplier=1024 * 5, base=1024 * 75
)
get_20gb_and_10gb_per_attempt = functools.partial(
    get_resources_per_gb, multiplier=1024 * 10, base=1024 * 20
)
get_30gb_and_10gb_per_attempt = functools.partial(
    get_resources_per_gb, multiplier=1024 * 10, base=1024 * 30
)
