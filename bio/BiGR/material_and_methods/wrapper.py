#!/usr/bin/python3.8
# -*- coding: utf-8 -*-

"""
Given a known list of parameters, this wrapper writes a mat&met for
BiGR pipelines
"""


from pathlib import Path


def get_version(path: Path, prefix: str, separator: str = "==") -> str:
    """Get tool version from environment path"""
    with path.open("r", enconding="utf-8") as env:
        for line in env:
            if line.startswith(prefix):
                return line[:-1].split(separator)[-1]


def fastp(params: str = "", wrappers_path: str = "") -> str:
    """Return Fastp Mat&Met"""
    environment_yaml_path = Path(
        wrappers_path, "bio", "fastp", "environment.yaml"
    )
    version = get_version(environment_yaml_path, "  - fastp", " =")
    mat_met = f"Fastq file cleaning has been performed by Fastp (version {version}), "

    if params == "":
        mat_met += "with default parameters. "
    else:
        mat_met += f"with the following optional parameters: '{params}'. "

    cite = 'Chen, Shifu, et al. "fastp: an ultra-fast all-in-one FASTQ preprocessor." Bioinformatics 34.17 (2018): i884-i890.'

    return mat_met, cite


def salmon(index: str = "",
           quant: str = "",
           wrappers_path: str = "",
           decoy: bool = True,
           in_house: bool = True) -> str:
    """Return Salmon Mat&Met"""
    environment_yaml_path = Path(
        wrappers_path, "bio", "salmon", "quant", "environment.yaml"
    )
    version = get_version(environment_yaml_path, "  - salmon")

    mat_met = f"Salmon (version {version}) was used for quantification step. "
    if decoy is True:
        mat_met += "Genome sequences were used as 'decoy sequences' while building transcriptome index. "

    if index == "":
        mat_met += "Default parameters were kept during indexation step. "
    else:
        mat_met += f"Salmon indexation step was done with the following optional parameters '{index}'. "

    if quant == "":
        mat_met += "Default parameters were kept during quantification step. "
    else:
        mat_met += f"Salmon quantification was done with the following optional parameters '{quant}'. "

    if in_house is True:
        mat_met += "In-house python scripts were used to aggregate and annotate single-sample counts. "

    cite = "Patro, Rob, et al. \"Salmon provides fast and bias-aware quantification of transcript expression.\" Nature methods 14.4 (2017): 417-419."

    return mat_met, cite

def snakemake(wrappers_path: str = "") -> str:
    """Return Snakemake Mat&Met"""
    environment_yaml_path = Path(
        wrappers_path, "bigr_pipelines", "common", "snakemake.yaml"
    )
    version = get_version(environment_yaml_path, "  - bioconda::snakemake")
    mat_met = f"This work was powered by Snakemake, version {version}."
    cite = "Köster, Johannes, and Sven Rahmann. \"Snakemake—a scalable bioinformatics workflow engine.\" Bioinformatics 28.19 (2012): 2520-2522."

    return mat_met, cite

def multiqc(wrappers_path: str = "") -> str:
    """Return MultiQC Mat&Met"""
    environment_yaml_path = Path(
        wrappers_path, "bio", "multiqc", "environment.yaml"
    )
    version = get_version(environment_yaml_path, "  - multiqc")
    mat_met = f"MultiQC (version {version}) was used to aggregate quality reports. "
    cite = "Ewels, Philip, et al. \"MultiQC: summarize analysis results for multiple tools and samples in a single report.\" Bioinformatics 32.19 (2016): 3047-3048."
    return mat_met, cite


def genome(organism: str, release: str, source: str, annotation: bool = True, variants: bool = False) -> str:
    """Cite organism"""
    mat_met = f" of {organism}, version {release}, from {source}, were used for this whole analysis. "
    if annotation and variants:
        mat_met = "The genome seuqence, annotation and known variants" + mat_met
    elif annotation:
        mat_met = "The genome seuqence, and annotation" + mat_met
    elif variants:
        mat_met = "The genome seuqence, and known variants" + mat_met
    else:
        mat_met = "The genome seuqences" + mat_met

    return mat_met


def bigr_salmon_quant_pipeline(organism: str = "Homo Sapiens",
                               release: str = "GRCh38.99",
                               source: str = "Ensembl",
                               fastp_params: str = "",
                               salmon_index_params: str = "",
                               salmon_quant_params: str = "",
                               decoy: bool = True,
                               annotation: bool = True):
    """Return Salmon quant pipeline Mat&Met"""
    snakemake_wrappers = "/mnt/beegfs/pipelines/snakemake-wrappers"
    mat_met = genome(
        organism=organism,
        release=release,
        source=source,
        annotation=annotation
    )
    cite = []

    m, c = fastp(params=fastp_params, wrappers_path=snakemake_wrappers)
    mat_met += m
    cite.appen(c)

    m, c = salmon(index=salmon_index_params, quant=salmon_quant_params, decoy=decoy, in_house=True, wrappers_path=snakemake_wrappers)
    mat_met += m
    cite.append(c)

    m, c = multiqc(wrappers_path=snakemake_wrappers)
    mat_met += m
    cite.append(c)

    m, c = snakemake(wrappers_path=snakemake_wrappers)
