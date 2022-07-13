#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Snakemake wrapper for SigProfiler Single Sample plots"""


import matplotlib.pyplot
import numpy
import pandas
import seaborn

from snakemake.utils import makedirs


from pathlib import Path
from scipy.spatial.distance import squareform
from scipy.spatial.distance import pdist, jaccard


# """
# Search, load, and transform data
# """
def read_sig_activities(path: Path) -> pandas.DataFrame:
    """Load SigProfiler results for a single sample"""
    print(f"Loading: {path}")
    return pandas.read_csv(path, sep="\t", header=0, index_col=0)


def load_all_sig_activities(*paths: Path, sort: bool = True) -> pandas.DataFrame:
    """Load multiple SigProfiler results and merge them"""
    print("Loading all signature activities")
    merged = None
    for path in paths:
        tmp = read_sig_activities(path)
        try:
            merged = pandas.merge(
                left = merged,
                right = tmp,
                how = "outer",
                left_index = True,
                right_index = True,
                sort = sort # Easier to read once sorted
                            # May be long for large tables/merges
            )
        except TypeError:
            # On the very first merge, it will fail due to "merged" being None
            merged = tmp

    merged.fillna(value=0, inplace=True)
    merged.index = [i.replace("Signature Subs-", "SBS") for i in merged.index]

    print("Activities merges:")
    print(merged.head())

    return merged


def percent(dataframe: pandas.DataFrame) -> pandas.DataFrame:
    """Convert column numbers in percentage"""
    for sample in dataframe.columns:
        total = dataframe[sample].sum()
        dataframe[sample] = [(i / total) * 100 for i in dataframe[sample]]

    return dataframe


def load_mutation_probabilities(path: str) -> pandas.DataFrame:
    """Load a mutation probability and format the table"""
    print(f"Loading probabilities for {path}")
    tmp = pandas.read_csv(path, sep="\t", header=0, index_col=None)
    sample = tmp["Sample Names"].tolist()[0]
    del tmp["Sample Names"]
    tmp.columns = [i.replace("Signature Subs-", "SBS") for i in tmp.columns]
    print(tmp.head())
    return tmp, sample


# """
# Build and display clustered heatmaps
# """
def clustermap_signatures(signatures: pandas.DataFrame,
                          png: str,
                          distance: str = "corr") -> None:
    """Save a clustered heatmap on disk"""
    print("Plotting clustermap")
    # Define color map
    cmap = seaborn.diverging_palette(
        h_neg=240,
        h_pos=10,
        as_cmap=True
    )

    # Build seaborn plot
    ax = None
    if distance == "corr":
        ax = seaborn.clustermap(
            signatures.corr(),
            cmap=cmap,
            method="average",
            metric="euclidean",
            robust=True,
            #annot=True
        )
    elif distance == "jaccard":
        ax = seaborn.clustermap(
            signatures.astype(bool),
            cmap=cmap,
            method="single",
            metric="jaccard",
            robust=True,
            #annot=True
        )

    # Rotate sample id to make them readable
    matplotlib.pyplot.setp(
        ax.ax_heatmap.yaxis.get_majorticklabels(),
        rotation=snakemake.params.get("ylabel_rotation", 0)
    )

    matplotlib.pyplot.setp(
        ax.ax_heatmap.xaxis.get_majorticklabels(),
        rotation=snakemake.params.get("xlabel_rotation", 90)
    )

    # Save result
    matplotlib.pyplot.savefig(png, bbox_inches="tight")
    matplotlib.pyplot.cla()
    matplotlib.pyplot.clf()
    matplotlib.pyplot.close()
    print("Done")


# """
# Build and display histograms
# """
def signatures_histogram(signatures: pandas.DataFrame, prefix: str) -> None:
    """Plot a histogram for a sample and its signatures"""
    print("Plotting all histograms")
    for sample_name in signatures.columns:
        print(f"Working on {sample_name}")
        # Define output png name
        png = ".".join([prefix, sample_name, "png"])
        if os.path.exists(png):
            print(f"{png} already exists, skipping ...")
            continue

        # Filter null values within a sample
        sample = signatures[signatures[sample_name].astype(bool)][sample_name]
        sig_name = sample.index.tolist()

        # Draw plot
        y_pos = numpy.arange(len(sample))
        fig, ax = matplotlib.pyplot.subplots()
        hbars = ax.barh(y_pos, sample)
        ax.set_yticks(y_pos, labels=sig_name)
        ax.invert_yaxis()  # labels read top-to-bottom
        ax.set_xlabel('Contribution')
        ax.set_title(f'Signatures in {sample_name}')

        # Save result
        matplotlib.pyplot.savefig(png, bbox_inches="tight")
        matplotlib.pyplot.cla()
        matplotlib.pyplot.clf()
        matplotlib.pyplot.close()
        print(f"Done for {sample_name}")
    print("Histograms done.")


def signatures_thistogram(signatures: pandas.DataFrame, prefix: str) -> None:
    """Plot a histogram for a sample and its signatures"""
    print("Plotting all transposed histograms")
    signatures = signatures.transpose()
    for signature in signatures.columns:
        print(f"Working on {signature}")
        # Define output png name
        png = ".".join([prefix, signature, "png"])
        if os.path.exists(png):
            print(f"{png} already exists, skipping ...")
            continue

        # Filter null values within a sample
        sig = signatures[signatures[signature].astype(bool)][signature]
        samples_name = sig.index.tolist()

        # Draw plot
        y_pos = numpy.arange(len(sig))
        fig, ax = matplotlib.pyplot.subplots()
        hbars = ax.barh(y_pos, sig)
        ax.set_yticks(y_pos, labels=samples_name)
        ax.invert_yaxis()  # labels read top-to-bottom
        ax.set_xlabel('Contribution')
        ax.set_title(f'Signatures in {signature}')

        # Save result
        matplotlib.pyplot.savefig(png, bbox_inches="tight")
        matplotlib.pyplot.cla()
        matplotlib.pyplot.clf()
        matplotlib.pyplot.close()
        print(f"Done for {signature}")
    print("Histograms done.")


# """
# Build and save Signature plots
# """
def plot_signatures_probabilities(prob: pandas.DataFrame, prefix: str) -> None:
    """Plot signature probabilities as an histogram"""
    for signature in prob.columns:
        if signature == "MutationTypes":
            continue

        png = ".".join([prefix, signature, "png"])
        if os.path.exists(png):
            print(f"{png} already exists, skipping ...")
            continue

        print(f"Working on {signature}'s probability -> {png}")

        # Build graph
        seaborn.set_theme(style="whitegrid")
        seaborn.set(rc={'figure.figsize':(30,10)})
        ax = seaborn.barplot(x="MutationTypes", y=signature, data=prob)
        ax.set_xticklabels(
            ax.get_xticklabels(),
            rotation=90,
            horizontalalignment='right'
        )

        # Save result
        matplotlib.pyplot.savefig(png, bbox_inches="tight", dpi=300)
        matplotlib.pyplot.cla()
        matplotlib.pyplot.clf()
        matplotlib.pyplot.close()
        print(f"Done for {signature}")

# """
# Main program and IO handlers
# """
signatures = None

if "matrix" in snakemake.output.keys():
    signatures = load_all_sig_activities(*snakemake.input["sigpro"])
    signatures.to_csv(
        snakemake.output["matrix"],
        sep="\t",
        header=True,
        index=True
    )

if "clustermap" in snakemake.output.keys():
    clustermap_signatures(signatures, snakemake.output["clustermap"])

if "jaccard_clustermap" in snakemake.output.keys():
    clustermap_signatures(
        signatures, snakemake.output["jaccard_clustermap"], distance="jaccard"
    )

if "histogram" in snakemake.output.keys():

    signatures_histogram(
        signatures, prefix=snakemake.params['hist_prefix'] + "absolute"
    )

if "transposed_hist" in snakemake.output.keys():
    signatures_thistogram(
        signatures,
        prefix=snakemake.params["hist_prefix"] + "absolute"
    )

signatures = percent(signatures) if isinstance(signatures, pandas.DataFrame) else None

if "percent_matrix" in snakemake.output.keys():
    signatures.to_csv(
        snakemake.output["percent_matrix"],
        sep="\t",
        header=True,
        index=True
    )

if "histogram" in snakemake.output.keys():
    signatures_histogram(
        signatures, prefix=snakemake.params['hist_prefix'] + "percent"
    )

if "transposed_hist" in snakemake.output.keys():
    makedirs(snakemake.output["transposed_hist"])
    signatures_thistogram(
        signatures,
        prefix=snakemake.params["hist_prefix"] + "percent"
    )


if "mut_prob" in snakemake.output.keys():
    makedirs(snakemake.output["mut_prob"])
    for prob_path in snakemake.input["prob"]:
        signatures, sample = load_mutation_probabilities(prob_path)
        prefix = snakemake.params["prob_prefix"] + sample
        plot_signatures_probabilities(signatures, prefix)
