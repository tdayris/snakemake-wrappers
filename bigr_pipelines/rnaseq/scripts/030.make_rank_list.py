#!/usr/bin/env python3
# coding: utf-8

import numpy
import logging
import pandas


def read_deseq2(
    path: str, padj_threshold: float = 1.0, fc_threshold: float = 0.0
) -> pandas.DataFrame:
    """Load DESeq2 dataset in memory"""
    logging.info("Loading %s", path)
    deseq = pandas.read_csv(path, sep="\t", header=0, index_col=0)
    deseq = deseq[deseq["padj"] <= padj_threshold]
    deseq = deseq[deseq["lfcSE"] >= fc_threshold]
    logging.debug(deseq.head())
    return deseq


def merge_tables(
    paths: List[str],
    on: str = "padj",
    padj_threshold: float = 1.0,
    fc_threshold: float = 0.0,
) -> pandas.DataFrame:
    """Merge n deseq2 tables"""
    path_iter = iter(paths)

    # Get the first data frame and convert padj in -log10(padj) if needed
    tmp = next(path_iter, None)
    tmp = read_deseq2(
        path=tmp, padj_threshold=padj_threshold, fc_threshold=fc_threshold
    )

    if on == "padj":
        tmp[on] = -numpy.log10(tmp[on])
    elif on == "log2FoldChange":
        on = "lfcSE"

    tmp = tmp[[on]]

    # Merge data frames one by one until the last one
    df = next(path_iter, None)
    while df is not None:
        df = read_deseq2(
            path=df, padj_threshold=padj_threshold, fc_threshold=fc_threshold
        )
        if on == "padj":
            df[on] = -numpy.log10(df[on])

        df = df[[on]]

        tmp = pandas.merge(
            left=tmp, right=df, left_index=True, right_index=True, how="outer"
        )
        df = next(path_iter, None)

    fill_value = 1.0
    if on == "padj":
        fill_value = 0.0
    tmp.fillna(fill_value)

    tmp.index.name = "ENSEMBL"

    return tmp


merged = merge_tables(
    paths=snakemake.input["deseq"],
    on=snakemake.params.get("rank_on", "padj"),
    padj_threshold=snakemake.params.get("padj_threshold", 0.05),
    fc_threshold=snakemake.params.get("fc_threshold", 0.7),
)

merged.to_csv(snakemake.output["tsv"], sep="\t", header=True, index=True)
