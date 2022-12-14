#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""Snakemake script to build multiqc additional figures from Salmon"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2019, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

import json
import logging

from typing import Dict, List, Union

Numbers = Union[int, float]

def get_sample_name(salmon_log: str) -> str:
    """
    Return sample name from log file path
    """
    return salmon_log.split("/")[-2]


def read_salmon_log(salmon_log: str, salmon_meta: str) -> Dict[str, Numbers]:
    """
    Read a salmon log file and builds output dictionary with
    number/percent of mapped fragments
    """
    result = {}

    logging.info(f"Reading {salmon_log}")
    with open(salmon_log, "r") as salmon_stream:
        data = json.load(salmon_stream)
    
    result["assigned_fragments"] = {
        "MSF": int(data["MSF"]),
        "OSF": int(data["OSF"]),
        "ISF": int(data["ISF"]),
        "MSR": int(data["MSR"]),
        "OSR": int(data["OSR"]),
        "ISR": int(data["ISR"]),
        "SF": int(data["SF"]),
        "SR": int(data["SR"]),
        "MU": int(data["MU"]),
        "OU": int(data["OU"]),
        "IU": int(data["IU"]),
        "U": int(data["U"]),
    }

    result["consistency"] = {
        "Consistent_Fragments": data["num_frags_with_concordant_consistent_mappings"],
        "Inconsistent_Fragments": data["num_frags_with_inconsistent_or_orphan_mappings"],
    }


    with open(salmon_meta, "r") as salmon_stream:
        data = json.load(salmon_stream)

    result["fragments_length"] = {
        "Length": data["frag_length_mean"],
        "SD": data["frag_length_sd"],
    }

    result["Percent_Mapped"] = data["percent_mapped"]

    result["fragments_num"] = {
        "Mapped": data["num_mapped"],
        "Decoy": data["num_decoy_fragments"],
        "Dovetail": data["num_dovetail_fragments"],
        "Filtered_VM": data["num_fragments_filtered_vm"],
        "Alignment_below_threshold": data["num_alignments_below_threshold_for_mapped_fragments_vm"],
    }

    return result

def main(
    salmon_logs: List[str], 
    salmon_metas: List[str], 
    salmon_mapping_plot: str,
    salmon_assigned_fragments: str,
    salmon_mapping_rates_general_stat: str
) -> None:
    """
    Read all salmon logs, save mapping rates in a multiqc config file
    """
    logging.info("Loading data...")
    data = {
        get_sample_name(salmon_log): read_salmon_log(salmon_log, salmon_meta)
        for salmon_log, salmon_meta in zip(salmon_logs, salmon_metas)
    }
    logging.info("Building bargraph config...")
    mqc_fragplot_config = {
        "id": "salmon_frag_type_count",
        "section_name": "Salmon Fragment assignation status",
        "description": "This histogram shows per-sample fragment assignation",
        "plot_type": "bargraph",
        "pconfig": {
            "id": "salmon_mapping_plot",
            "title": "Salmon fragment assignation status",
            "xlab": "Samples",
            "ylab": "Number of fragments",
        },
        "data": {
            sample_name: sample_data["assigned_fragments"]
            for sample_name, sample_data in data.items()
        }
    }

    with open(salmon_mapping_plot, "w") as salmon_mapping_plot_stream:
        json.dump(obj=mqc_fragplot_config, fp=salmon_mapping_plot_stream)
    logging.info(f"Salmon fragment assignation qc build: {salmon_mapping_plot}")

    mqc_assignplot_config = {
        "id": "salmon_assigned_fragments",
        "section_name": "Salmon Fragment assignation type",
        "description": "This histogram shows per-sample fragment mapping counts",
        "plot_type": "bargraph",
        "pconfig": {
            "id": "salmon_mapping_plot",
            "title": "Salmon fragment assignation counts",
            "xlab": "Samples",
            "ylab": "Number of fragments",
        },
        "data": {
            sample_name: sample_data["fragments_num"]
            for sample_name, sample_data in data.items()
        }
    }

    with open(salmon_assigned_fragments, "w") as salmon_mapping_plot_stream:
        json.dump(obj=mqc_assignplot_config, fp=salmon_mapping_plot_stream)
    logging.info(f"Salmon fragment count qc build: {salmon_assigned_fragments}")

    logging.info("Building general stats config...")
    mqc_general_stat = {
        "id": "salmon_mapping_rates_general_stat",
        "plot_type": "generalstats",
        "pconfig": {
            "Percent_Mapped": {
                "max": 100,
                "min": 0,
                "suffix": "%",
            },
        },
        "data": {
            sample_name: sample_data["Percent_Mapped"]
            for sample_name, sample_data in data.items()
        }
    }

    with open(salmon_mapping_rates_general_stat, "w") as salmon_mapping_plot_stream:
        json.dump(obj=mqc_general_stat, fp=salmon_mapping_plot_stream)
    logging.info(f"Salmon general table qc build: {salmon_mapping_rates_general_stat}")


if __name__ == "__main__":
    logging.basicConfig(
        filename=snakemake.log[0],
        filemode="w", 
        level=logging.DEBUG
    )

    try:
        main(
            salmon_logs=snakemake.input["salmon_logs"],
            salmon_metas=snakemake.input["salmon_meta"],
            salmon_mapping_plot=snakemake.output["mapping_status_mqc"],
            salmon_assigned_fragments=snakemake.output["mapping_counts_mqc"],
            salmon_mapping_rates_general_stat=snakemake.output["salmon_general_table"],
        )
    except Exception as e:
        logging.error(e)
        raise