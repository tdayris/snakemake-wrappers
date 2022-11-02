# /usr/bin/env python3
# coding: utf-8

import logging
from operator import index
import os.path
import pandas
import seaborn
import matplotlib.pyplot

from functools import partial
from pathlib import Path


from typing import Any, Callable, Dict, Generator, Optional, Union

PathType = Union[str, Path]


def search_logs(path: PathType) -> Generator[Dict[str, Union[Path, int]], None, None]:
    """
    Search for valid logging files
    (aka. with /usr/bin/time content)
    """
    logging.info("Searching for /usr/bin/time reports in %s", str(path.absolute()))
    if isinstance(path, str):
        path = Path(path)

    for file in path.iterdir():
        line_number = 0

        if not file.name.endswith(".err"):
            continue

        with file.open("r") as error_log:
            line_iter = iter(error_log)
            line = next(line_iter, None)
            line_number += 1
            while line is not None:
                last_line = line
                line = next(line_iter, None)

        if last_line.startswith("	Exit status: "):
            yield {"path": file, "start": line_number - 22, "end": line_number}


def get_value(
    line: str,
    separator: str = ": ",
    to: Any = float,
    strip: Optional[str] = None,
) -> Any:
    """Get the value of a line"""
    raw_value = line.split(sep=separator)[-1]
    if strip is not None:
        raw_value = raw_value.strip(strip)
    return to(raw_value)


split_percent = partial(get_value, strip="%")
split_float = partial(get_value)
split_string = partial(get_value, to=str)


def parse_error_log(path: str, start: int, end: int) -> Dict[str, Union[float, str]]:
    """Load error log information of interest in memory"""
    logging.info("Parsing: %s", path)
    content = {
        "date": "Unknown",
        "Reserved_threads": 1,
        "Reserved_memory": 1,
        "CPU": 0.0,
        "time": "0:0.0",
        "max_memory": 1.0,
        "name": "Unknown",
        "jobid": 0,
        "node": "Unknown",
        "step": 0,
    }

    with open(path, "r") as error_log:
        for idx, line in enumerate(error_log):
            line = line[:-1]
            if line.startswith("Provided cores: 2") and idx < 16:
                content["Reserved_threads"] = split_float(line)
            elif line.startswith("[") and line.endswith("]") and idx < 16:
                content["date"] = line.strip("[]")
            elif line.startswith("    threads: "):
                content["Reserved_threads"] = split_float(line)
            elif line.startswith("    resources: mem_mb"):
                content["Reserved_memory"] = float(
                    line.split(": ")[1].split(", ")[0].split("=")[1]
                )
            elif line.startswith("	Percent of CPU this job got:"):
                content["CPU"] = split_percent(line)
            elif line.startswith("	Elapsed (wall clock) time (h:mm:ss or m:ss):"):
                content["time"] = split_string(line)
            elif line.startswith("	Maximum resident set size (kbytes): "):
                content["max_memory"] = (
                    split_float(line) / 1_024
                )  # Convert kbytes to mbytes

    chomp = str(path.name).split(".")
    content["name"] = chomp[0].split("-")[-1]
    content["jobid"] = int(chomp[4].split("-")[1])
    content["node"] = chomp[4].split("-")[-1]
    content["step"] = int(chomp[1])

    content["Memory_Usage"] = content["max_memory"] / content["Reserved_memory"]
    content["CPU_Usage"] = (content["CPU"] / 100) * content["Reserved_threads"]

    return content


logging.basicConfig(level=logging.DEBUG)
parsed_logs_path = "logs/parsed.logs.tsv"
barplot = True

if not os.path.exists(parsed_logs_path):
    logging.info("Creating the %s file", parsed_logs_path)
    jobs_dict = {}

    for log in search_logs(Path("logs/slurm/")):
        err = parse_error_log(log["path"], log["start"], log["end"])
        jobs_dict[f"{err['name']}.{err['jobid']}"] = err

    jobs = pandas.DataFrame.from_dict(data=jobs_dict,orient="index",)
    
    jobs.index.name = "Job_name"
    jobs.to_csv(parsed_logs_path, sep="\t", header=True, index=True)
else:
    logging.info("Parsed logs found: %s", parsed_logs_path)
    jobs = pandas.read_csv(parsed_logs_path, sep="\t", header=0, index_col=0)

logging.debug(jobs.head())

if barplot is True:
    png_out = "MemoryConsumption.png"
    seaborn.catplot(data=jobs, x="name", y="max_memory", kind="bar")
    matplotlib.pyplot.xticks(rotation=90)

    logging.info("Gene plot saved to %s", png_out)
    matplotlib.pyplot.savefig(
        png_out,
        bbox_inches="tight"
    )