#!/usr/bin/python3
# -*- coding: utf-8 -*-

import json
import logging

from typing import Dict, List

def get_sample_name(salmon_log: str) -> str:
    """
    Return sample name from log file path
    """
    return salmon_log.split("/")[-2]


def read_salmon_log(salmon_log: str, salmon_meta: str) -> Dict[str, int]:
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
        "Inconsistent_Fregments": data["num_frags_with_inconsistent_or_orphan_mappings"],
    }


    with open(salmon_meta, "r") as salmon_stream:
        data = json.load(salmon_stream)

    result["fragments_length"] = {
        "Length": data["frag_length_mean"],
        "SD": data["frag_length_sd"],
    }

    result["percent"] = {
        "Percent_Mapped": data["percent_mapped"],
    }

    result["fragments_num"] = {
        "Mapped": data["num_mapped"],
        "Decoy": data["num_decoy_fragments"],
        "Dovetail": data["num_dovetail_fragments"],
        "Filtered_VM": data["num_fragments_filtered_vm"],
        "Alignment_below_threshold": data["num_alignments_below_threshold_for_mapped_fragments_vm"],
    }

    return result

def main(salmon_logs: List[str], salmon_metas: List[str], outfile: str) -> None:
    """
    Read all salmon logs, save mapping rates in a multiqc config file
    """
    logging.info("Loading data...")
    data = {
        get_sample_name(salmon_log): read_salmon_log(salmon_log, salmon_meta)
        for salmon_log, salmon_meta in zip(salmon_logs, salmon_metas)
    }
    logging.info("Building bargraph config...")
    mqc_plot_config = {
        "id": "salmon_mapping_rates",
        "section_name": "Salmon Mapping Rates",
        "description": "This histogram shows per-sample mapping rates",
        "plot_type": "bargraph",
        "pconfig": {
            "id": "salmon_mapping_plot",
            "title": "Salmon mapping rates",
            "xlab": "Samples",
            "ylab": "Number of reads",
        },
        "data": {
            sample_name: sample_data["assigned_fragments"]
            for sample_name, sample_data in data.items()
        }
    }

    logging.info("Building general stats config...")
    mqc_general_stat = {
        "id": "salmon_mapping_rates_general_stat",
        "plot_type": "generalstats",
        "pconfig": {
            "Mapping_Rate": {
                "max": 100,
                "min": 0,
                "suffix": "%",
            },
            "Consistent_Mappings": {
                "max": 100,
                "min": 0,
                "suffix": "%",
            },
            "Inconsistent_Mappings": {
                "max": 100,
                "min": 0,
                "suffix": "%",
            },
        },
        "data": {
            sample_name: sample_data["percent"]
            for sample_name, sample_data in data.items()
        }
    }


if __name__ == "__main__":
    logging.basicConfig(
        filename=snakemake.log[0],
        filemode="w", 
        level=logging.DEBUG
    )