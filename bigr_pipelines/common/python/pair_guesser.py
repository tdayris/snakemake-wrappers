#!/usr/bin/env python3
# coding: utf-8

"""Read paired dataset paths from odin kdi procedures, and guess splitted samples"""

import argparse
import logging
import re
import sys

from typing import Any, Dict, List, Optional


class CustomFormatter(logging.Formatter):

    cyan = "\033[1;36m"
    orange = "\x1b[33;20m"
    red = "\x1b[31;20m"
    green = "\x1b[31;1m"
    reset = "\x1b[0m"
    format = "%(levelname)s - %(asctime)s - %(name)s - %(message)s (%(filename)s:%(lineno)d)"

    FORMATS = {
        logging.DEBUG: green + format + reset,
        logging.INFO: cyan + format + reset,
        logging.WARNING: orange + format + reset,
        logging.ERROR: red + format + reset,
        logging.CRITICAL: red + format + reset
    }

    def format(self, record):
        log_fmt = self.FORMATS.get(record.levelno)
        formatter = logging.Formatter(log_fmt)
        return formatter.format(record)



def cmd_parser(args: Any) -> argparse.Namespace:
    """Parse command line arguments"""
    logging.info("Parsing command line...")
    parser = argparse.ArgumentParser(
        description="",
        epilog="This script does not perform any magic. Check results."
    )

    parser.add_argument(
        "-i", "--input",
        help="Path to input TSV file",
        default="paired_dataset_paths.txt",
        type=str
    )

    parser.add_argument(
        "-d", "--design",
        help="Path to the output design file with merged fastq files",
        default="design.tsv",
        type=str
    )

    parser.add_argument(
        "-e", "--extension",
        help="Space separated list of suffixes to remove to guess sample name",
        default=[".gz", ".fq", ".fastq", "_001", "_R1"],
        nargs="+"
    )

    return parser.parse_args(args)


def try_remove_ext(filename: str, possible_exts: Optional[List[str]] = None) -> str:
    """From  a list of possible extensions, remove them one after each other"""
    if possible_exts is None:
        possible_exts = [
            ".gz",
            ".fastq",
            ".fq",
            "_001",
            "_R1",
        ]

    for ext in possible_exts:
        filename = re.sub(ext, "", filename)
    return filename


def index_paired_dataset(path: str, possible_exts: Optional[List[str]]) -> Dict[str, Dict[str, List[str]]]:
    logging.info("Indexing input file...")
    result = {}
    with open(path, "r") as paired_dataset:
        for line in paired_dataset:
            left, right = line[:-1].split("\t")
            sample_name = try_remove_ext(left.split("/")[-1], possible_exts)
            try:
                result[sample_name]["R1"].append(left)
                result[sample_name]["R2"].append(right)
            except KeyError:
                result[sample_name] = {"R1": [left], "R2": [right]}

    return result


def create_design_file(path: str, index: Dict[str, Dict[str, List[str]]]) -> None:
    logging.info("Saving design file...")
    with open(path, "w") as design:
        design.write("Sample_id\tUpstream_file\tDownstream_file\n")
        for sample in index.keys():
            if len(index[sample]["R1"]) > 1:
                logging.info("{} seems to be divided in {} files".format(sample, len(index[sample]["R1"])))
            line = "\t".join([
                sample,
                ",".join(index[sample]["R1"]),
                ",".join(index[sample]["R2"]),
            ]) + "\n"
            design.write(line)


if __name__ == "__main__":
    # create logger with 'spam_application'
    logger = logging.getLogger("pair_guesser")
    logger.setLevel(logging.DEBUG)

    # create console handler with a higher log level
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)

    ch.setFormatter(CustomFormatter())

    logger.addHandler(ch)

    cmd = cmd_parser(sys.argv[1:])
    logging.debug(cmd)

    indexed_input = index_paired_dataset(cmd.input, cmd.extension)

    create_design_file(cmd.design, indexed_input)