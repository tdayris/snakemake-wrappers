#!/usr/bin/python3
# coding: utf-8

"""
A snakemake wrapper that builts a set of yaml files from a snakemake extended-benchmark set of directories.
This script expects the benchmark directory to be formetted as follows: root_benchmark_dir / rule_name / reports
"""

import ast
import argparse
import dataclasses
import humanfriendly
import numpy
import os
import pandas
import pathlib
import psutil
import rich_argparse
import sklearn.linear_model
import shutil
import sys
import typing
import yaml

from loguru import logger


@dataclasses.dataclass
class MockSizeEstimate:
    size_mb: int


def fmt_bytes(bt: int) -> str:
    return humanfriendly.format_size(bt if not numpy.isnan(bt) else 0)

def fmt_mb(mb: int) -> str:
    return fmt_bytes(mb * 10 ** 6)

def fmt_gb(gb: int) -> str:
    return fmt_mb(gb * 10 ** 3)

def fmt_seconds(sec: int) -> str:
    return humanfriendly.format_timespan(sec if not numpy.isnan(sec) else 0)

def fmt_minutes(minutes: int) -> str:
    return fmt_seconds(minutes * 60)

def fmt_hours(hours: int) -> str:
    return fmt_minutes(hours * 60)

def fmt_days(days: int) -> str:
    return fmt_hours(days * 24)

def fmt_month(month: int) -> str:
    return fmt_days(month * 30)

def compute_ratio(num: int, den: int) -> float:
    if any(v == 0 for v in [num, den]):
        return numpy.nan
    return num / den

def series_of_dicts_to_frame(series: pandas.Series, series_name: str) -> pandas.DataFrame:
    """Expand a pandas column with dict/json into a pandas DataFrame"""
    records: pandas.DataFrame = pandas.json_normalize(series.dropna().map(ast.literal_eval))
    records.columns = [f"{series_name}_{c}" for c in records.columns]
    return records

def train_linear_regression(
    io_in: pandas.DataFrame,
    target: pandas.DataFrame,
) -> dict[str, sklearn.linear_model.LinearRegression | float]:
    """Build model to predict best resources based on what's been seen before"""
    model = sklearn.linear_model.LinearRegression().fit(io_in, target)
    return {
        "model": model, 
        "coef": float(model.coef_[0][0]),
        "intercept": float(model.intercept_[0]),
    }

def predict_resources(
    model: sklearn.linear_model.LinearRegression, 
    io_in: float,
) -> dict[str, float]:
    """Predict rss, seconds, out from a trained model"""
    pred: numpy.ndarray[tuple[int, ...], numpy.dtype[numpy.float64]] = numpy.exp(
        model.predict([[numpy.log(io_in)]])[0]
    )
    return dict(zip(["max_rss", "seconds", "io_out"], pred))


def process_extended_benchmark(
    data: pandas.DataFrame,
) -> pandas.DataFrame:
    """
    Load a benchmark table in memory, filter and builds required columns
    """
    logger.info("Enhancing benchmark data")
    benchmark: pandas.DataFrame = data.copy()
    cols_to_keep: list[str] = [
        "s",
        "max_rss",
        "io_in",
        "io_out",
        "jobid",
        "rule_name",
        "threads",
        "cpu_usage",
        "resources",
    ]

    # Expand resources and evaluate efficiency
    resources: pandas.DataFrame = series_of_dicts_to_frame(
        series=benchmark.resources.copy(),
        series_name="resources",
    )
    benchmark = pandas.concat([benchmark, resources], axis=1)
    benchmark["memory_efficiency"] = benchmark["max_rss"] / benchmark["resources_mem_mb"]
    benchmark["time_efficiency"] = benchmark["s"] / (benchmark["resources_runtime"] * 60)
    benchmark["disk_efficiency"] = ( benchmark["io_out"] + benchmark["io_in"] ) / benchmark["resources_disk_mb"]

    # Expend input data size
    input_data_size: pandas.DataFrame = series_of_dicts_to_frame(
        series=benchmark.input_size_mb.copy(),
        series_name="input_size_mb",
    )
    input_data_size["input_size_mb_total"] = input_data_size.sum(axis=1)
    input_data_size["input_size_gb_total"] = input_data_size["input_size_mb_total"] / 1_000
    benchmark = pandas.concat([benchmark, input_data_size], axis=1)

    # Always add 10 additional minutes to have a small security margin
    benchmark["seconds"] = benchmark["s"]
    benchmark["minutes"] = numpy.ceil(benchmark["seconds"] / 60).astype(numpy.uint16)

    # Always add 300 additional Mb  to have a small security margin
    benchmark["memory_mb"] = benchmark["max_rss"]
    benchmark["memory_gb"] = [max(1, rss / 1_000) for rss in benchmark["memory_mb"]]

    # Always add 1gb to input size to have a small security margin
    benchmark["io_in"] = [max(1, i) for i in benchmark["io_in"]]
    benchmark["input_mb"] = benchmark["io_in"] / 1_000_000
    benchmark["input_gb"] = numpy.ceil(benchmark["input_mb"] / 1_000).astype(numpy.uint16)
    benchmark["io_out"] = [max(1, o) for o in benchmark["io_out"]]
    benchmark["output_mb"] = benchmark["io_out"] / 1_000_000
    benchmark["output_gb"] = numpy.ceil(benchmark["output_mb"] / 1_000).astype(numpy.uint16)
    benchmark["disk_bytes"] = benchmark["io_in"] + benchmark["io_out"]
    benchmark["disk_mb"] = benchmark["disk_bytes"] / 1_000_000
    benchmark["disk_gb"] = benchmark["disk_mb"] / 1_000

    # Convert to human readable metrics
    benchmark["human_seconds"] = [fmt_seconds(s) for s in benchmark["s"]]
    benchmark["human_memory"] = [fmt_mb(mb) for mb in benchmark["max_rss"]]
    benchmark["human_input"] = [fmt_bytes(bt) for bt in benchmark["io_in"]]
    benchmark["human_output"] = [fmt_bytes(bt) for bt in benchmark["io_out"]]

    # Compute ratios, encouting for null values
    benchmark["time_ratio"] = [compute_ratio(time, size) for time, size in zip(benchmark.seconds, benchmark.io_in)]
    benchmark["ram_ratio"] = [compute_ratio(memory_mb, size) for memory_mb, size in zip(benchmark.memory_mb, benchmark.io_in)]
    benchmark["out_ratio"] = [compute_ratio(out, size) for out, size in zip(benchmark.io_out, benchmark.io_in)]

    cols_to_keep += [
        "seconds",
        "minutes",
        "human_seconds",
        "memory_mb",
        "memory_gb",
        "human_memory",
        "input_mb",
        "input_gb",
        "human_input",
        "output_mb",
        "output_gb",
        "disk_bytes",
        "disk_mb",
        "disk_gb",
        "input_size_mb_total",
        "input_size_gb_total",
        "human_output",
        "time_ratio",
        "ram_ratio",
        "out_ratio",
        "memory_efficiency",
        "time_efficiency",
        "disk_efficiency",
    ]
    benchmark = benchmark[cols_to_keep]
    logger.debug(benchmark)
    return benchmark


def cat_dataframes(
    paths: list[str | pathlib.Path],
) -> pandas.DataFrame:
    """
    Load multiple TSV in memory and concatenate them
    """
    logger.info("Loading and concatenating benchmark files")
    frames = [
        pandas.read_csv(path, sep="\t", index_col=None, header=0) for path in paths
    ]
    return pandas.concat(frames, ignore_index=True)


def compute_delta(base_res: float | int, max_res: float | int, input_size: float | int, attempt: str) -> str:
    """Return best delta given input size and reservations"""
    if any(x == 0.0 for x in (base_res, max_res, input_size)):
        if base_res < max_res:
            return f"{attempt} * {max(max_res - base_res, 0.1):.2f} + {base_res:.2f}"
        return f"{attempt} * {base_res:.2f} + {base_res:.2f}"
    return f"{attempt} * {max(max_res / input_size, 0.1):.2f} * input.size_mb + {base_res:.2f}"


def resource_function(
    known: pandas.DataFrame,
    resource_type: str = "time",
    mode: str = "slurm",
) -> dict[str, str]:
    """Use a simple method to assess resrouces requirement"""
    min_reservation: float | int = 0.0
    base_reservation: float | int = 0.0
    max_reservation: float | int = 0.0
    delta: float | int = 0.0
    max_allowed_reservation: float | int = 365
    min_allowed_reservation: float | int = 0
    attempt: str = "attempt"
    input_size: float | int = known.input_size_mb_total.mean()

    if resource_type == "time":
        if mode == "slurm":
            min_reservation = known.minutes.min()
            base_reservation = known.minutes.quantile(0.75)
            input_size = known.input_size_mb_total.quantile(0.75)
            max_reservation = known.minutes.max()
            # Roughly 30 days in minutes
            max_allowed_reservation = 43800
            min_allowed_reservation = 10
        elif mode == "google":
            min_reservation = known.seconds.min()
            base_reservation = known.seconds.quantile(0.75)
            input_size = known.input_size_mb_total.quantile(0.75)
            max_reservation = known.seconds.max()
            # Roughly 30 days in seconds
            max_allowed_reservation = 2628000
            min_allowed_reservation = 600
        elif mode == "local":
            min_reservation = known.minutes.min()
            base_reservation = known.minutes.quantile(0.75)
            input_size = known.input_size_mb_total.quantile(0.75)
            max_reservation = known.minutes.max()
            # Roughly 300 days in minutes
            max_allowed_reservation = 438000
            min_allowed_reservation = 0

    elif resource_type == "ram":
        if mode == "slurm":
            min_reservation = known.memory_mb.min()
            base_reservation = known.memory_mb.quantile(0.75)
            input_size = known.input_size_mb_total.quantile(0.75)
            max_reservation = known.memory_mb.max()
            attempt = "(attempt  - 1)"
            # Max 150Gb
            max_allowed_reservation = 150000
            min_allowed_reservation = 750
        elif mode == "google":
            min_reservation = known.memory_gb.min()
            base_reservation = known.memory_gb.quantile(0.75)
            input_size = known.input_size_mb_total.quantile(0.75)
            max_reservation = known.memory_gb.max()
            attempt = "(attempt  - 1)"
            # Max 1To, min 3Go for VM
            max_allowed_reservation = 100000
            min_allowed_reservation = 3
        elif mode == "local":
            min_reservation = known.memory_mb.min()
            base_reservation = known.memory_mb.quantile(0.75)
            input_size = known.input_size_mb_total.quantile(0.75)
            max_reservation = known.memory_mb.max()
            attempt = "(attempt  - 1)"
            # Max set to max installed RAM, min is 1Mb
            max_allowed_reservation = psutil.virtual_memory().total
            min_allowed_reservation = 1

    elif resource_type == "disk":
        script_dir: pathlib.Path = pathlib.Path(__file__).resolve().parent
        if mode == "slurm":
            base_reservation = known.disk_mb.mean()
            max_reservation = known.disk_mb.max()
            min_reservation = known.disk_mb.min()
            input_size = known.input_size_mb_total.quantile(0.75)
            # Max is max available space on disk, min is 2Go.
            max_allowed_reservation = shutil.disk_usage(script_dir).free
            min_allowed_reservation = 2
        elif mode == "google":
            min_reservation = known.disk_gb.min()
            base_reservation = known.disk_gb.mean()
            max_reservation = known.disk_gb.max()
            input_size = known.disk_gb.quantile(0.75)
            # Max is 2To, min is 30 GB.
            max_allowed_reservation = 2000
            min_allowed_reservation = 30
        elif mode == "local":
            min_reservation = known.disk_mb.min()
            base_reservation = known.disk_mb.mean()
            max_reservation = known.disk_mb.max()
            input_size = known.input_size_mb_total.quantile(0.75)
            # Max is max available space on disk, min is 0.
            max_allowed_reservation = shutil.disk_usage(script_dir).free
            min_allowed_reservation = 0

    delta = compute_delta(base_reservation, max_reservation, input_size, attempt)
    anchored_function = f"round(min({max_allowed_reservation}, max({min_allowed_reservation}, {delta})))"

    resource_function: str = str(anchored_function)
    if resource_type == "time" and mode == "google":
        resource_function = f"str({anchored_function}) + 's'"
    estimate = {
        f"{resource_type}_function": resource_function,
    }

    input = MockSizeEstimate(size_mb=1_000)
    if resource_type == "time":
        if mode == "google":
            attempt: int = 1
            estimate[f"raw_{resource_type}_attempt1"] = eval(anchored_function)
            estimate[f"{resource_type}_attempt1"] = fmt_seconds(estimate[f"raw_{resource_type}_attempt1"])
            attempt: int = 2
            estimate[f"raw_{resource_type}_attempt2"] = eval(anchored_function)
            estimate[f"{resource_type}_attempt2"] = fmt_seconds(estimate[f"raw_{resource_type}_attempt2"])
            attempt: int = 13
            estimate[f"raw_{resource_type}_attempt13"] = eval(anchored_function)
            estimate[f"{resource_type}_attempt13"] = fmt_seconds(estimate[f"raw_{resource_type}_attempt13"])
        else:
            attempt: int = 1
            estimate[f"raw_{resource_type}_attempt1"] = eval(anchored_function)
            estimate[f"{resource_type}_attempt1"] = fmt_minutes(estimate[f"raw_{resource_type}_attempt1"])
            attempt: int = 2
            estimate[f"raw_{resource_type}_attempt2"] = eval(anchored_function)
            estimate[f"{resource_type}_attempt2"] = fmt_minutes(estimate[f"raw_{resource_type}_attempt2"])
            attempt: int = 13
            estimate[f"raw_{resource_type}_attempt13"] = eval(anchored_function)
            estimate[f"{resource_type}_attempt13"] = fmt_minutes(estimate[f"raw_{resource_type}_attempt13"])
    else:
        if mode == "google":
            attempt: int = 1
            estimate[f"raw_{resource_type}_attempt1"] = eval(anchored_function)
            estimate[f"{resource_type}_attempt1"] = fmt_gb(estimate[f"raw_{resource_type}_attempt1"])
            attempt: int = 2
            estimate[f"raw_{resource_type}_attempt2"] = eval(anchored_function)
            estimate[f"{resource_type}_attempt2"] = fmt_gb(estimate[f"raw_{resource_type}_attempt2"])
            attempt: int = 13
            estimate[f"raw_{resource_type}_attempt13"] = eval(anchored_function)
            estimate[f"{resource_type}_attempt13"] = fmt_gb(estimate[f"raw_{resource_type}_attempt13"])
        else:
            attempt: int = 1
            estimate[f"raw_{resource_type}_attempt1"] = eval(anchored_function)
            estimate[f"{resource_type}_attempt1"] = fmt_mb(estimate[f"raw_{resource_type}_attempt1"])
            attempt: int = 2
            estimate[f"raw_{resource_type}_attempt2"] = eval(anchored_function)
            estimate[f"{resource_type}_attempt2"] = fmt_mb(estimate[f"raw_{resource_type}_attempt2"])
            attempt: int = 13
            estimate[f"raw_{resource_type}_attempt13"] = eval(anchored_function)
            estimate[f"{resource_type}_attempt13"] = fmt_mb(estimate[f"raw_{resource_type}_attempt13"])

    return estimate


def rule_profile(
    known: pandas.DataFrame,
    machine_types_path: str | pathlib.Path,
    mode: str = "local",
    googlebatch_boot_disk_type: str = "hyperdisk-balanced",
) -> list[str]:
    logger.info(f"Building rule profile for {known['rule_name'].tolist()[0]}, {mode=}")
    data: dict[str, str | int] = {
        **resource_function(known, "time", mode),
        **resource_function(known, "ram", mode),
        **resource_function(known, "disk", mode),
    }
    data["rule_name"] = known["rule_name"].tolist()[0]
    text: list[str] = [f"  {data['rule_name']}:", f"    # Threads: {known['threads'].tolist()[0]}"]
    if mode == "local":
            
        text += [
            f"    # Max time seen for this rule: {fmt_seconds(known['s'].max())}",
            f"    # Mean time seen for this rule: {fmt_seconds(known['s'].mean())} ± {fmt_seconds(known['s'].std())}",
            f"    # Time reservation for 1st attempt: {data['time_attempt1']}",
            f"    # Time reservation for 2nd attempt: {data['time_attempt2']}",
            f"    # Time reservation for 13th attempt: {data['time_attempt13']}",
            f"    runtime: {data['time_function']}",
            f"    # Max ram seen for this rule {fmt_mb(known['max_rss'].max())}",
            f"    # Mean RAM seen for this rule: {fmt_mb(known['memory_mb'].mean())} ± {fmt_mb(known['memory_mb'].std())}",
            f"    # RAM reservation for 1st attempt: {data['ram_attempt1']}",
            f"    # RAM reservation for 2nd attempt: {data['ram_attempt2']}",
            f"    # RAM reservation for 13th attempt: {data['ram_attempt13']}",
            f"    mem: {data['ram_function']}",
            f"    # Max IO seen for this rule {fmt_mb(known['disk_bytes'].max())}",
            f"    # Mean disk used for this rule: {fmt_bytes(known['disk_bytes'].mean())} ± {fmt_bytes(known['disk_bytes'].std())}",
            f"    # Disk reservation for 1st attempt: {data['disk_attempt1']}",
            f"    # Disk reservation for 2nd attempt: {data['disk_attempt2']}",
            f"    # Disk reservation for 13th attempt: {data['disk_attempt13']}",
            f"    disk_mb: {data['disk_function']}",
        ]
    elif mode == "slurm":
        text += [
            f"    # Max time seen for this rule: {fmt_seconds(known['s'].max())}",
            f"    # Mean time seen for this rule: {fmt_seconds(known['s'].mean())} ± {fmt_seconds(known['s'].std())}",
            f"    # Time reservation for 1st attempt: {data['time_attempt1']}",
            f"    # Time reservation for 2nd attempt: {data['time_attempt2']}",
            f"    # Time reservation for 13th attempt: {data['time_attempt13']}",
            f"    runtime: {data['time_function']}",
            f"    # Max ram seen for this rule {fmt_mb(known['max_rss'].max())}",
            f"    # Mean RAM seen for this rule: {fmt_mb(known['memory_mb'].mean())} ± {fmt_mb(known['memory_mb'].std())}",
            f"    # RAM reservation for 1st attempt: {data['ram_attempt1']}",
            f"    # RAM reservation for 2nd attempt: {data['ram_attempt2']}",
            f"    # RAM reservation for 13th attempt: {data['ram_attempt13']}",
            f"    mem: {data['ram_function']}",
            f"    mem_mb: {data['ram_function']}",
            f"    # Max IO seen for this rule {fmt_bytes(known['disk_bytes'].max())}",
            f"    # Mean disk used for this rule: {fmt_bytes(known['disk_bytes'].mean())} ± {fmt_bytes(known['disk_bytes'].std())}",
            f"    # Disk reservation for 1st attempt: {data['disk_attempt1']}",
            f"    # Disk reservation for 2nd attempt: {data['disk_attempt2']}",
            f"    # Disk reservation for 13th attempt: {data['disk_attempt13']}",
            f"    disk_mb: {data['disk_function']}",
            f"    # Cluster account: sacctmgr -nP show assoc where user=$USER format=account",
            f"    slurm_account: 'igr_hpc'",
            f"    # Cluster name: scontrol show config | grep ClusterName",
            f"    cluster: 'flamingo'",
        ]
    elif mode == "google":
        data["googlebatch_machine_types"], data["googlebatch_machine_prices"] = chose_google_machine_type(
            known=known,
            machine_types_path=machine_types_path,
        )

        text += [
            f"    # Max time seen for this rule: {fmt_seconds(known['s'].max())}",
            f"    # Mean time seen for this rule: {fmt_seconds(known['s'].mean())} ± {fmt_seconds(known['s'].std())}",
            f"    # Time reservation for 1st attempt: {data['time_attempt1']}, about {(data['raw_time_attempt1'] / 60 / 60) / data['googlebatch_machine_prices'][0]:.2f} $",
            f"    # Time reservation for 2nd attempt: {data['time_attempt2']}, about {(data['raw_time_attempt2'] / 60 / 60) / data['googlebatch_machine_prices'][0]:.2f} $",
            f"    # Time reservation for 13th attempt: {data['time_attempt13']}, about {(data['raw_time_attempt13'] / 60 / 60) / data['googlebatch_machine_prices'][0]:.2f} $",
            f"    googlebatch_max_run_duration: {data['time_function']}",
            f"    # Max ram seen for this rule {fmt_mb(known['max_rss'].max())}",
            f"    # Mean RAM seen for this rule: {fmt_mb(known['memory_mb'].mean())} ± {fmt_mb(known['memory_mb'].std())}",
            f"    # RAM reservation for 1st attempt: {data['ram_attempt1']}",
            f"    # RAM reservation for 2nd attempt: {data['ram_attempt2']}",
            f"    # RAM reservation for 13th attempt: {data['ram_attempt13']}",
            f"    googlebatch_memory: {data['ram_function']}",
            f"    # Max IO seen for this rule {fmt_bytes(known['disk_bytes'].max())}",
            f"    # Mean disk used for this rule: {fmt_bytes(known['disk_bytes'].mean())} ± {fmt_bytes(known['disk_bytes'].std())}",
            f"    # Disk reservation for 1st attempt: {data['disk_attempt1']}",
            f"    # Disk reservation for 2nd attempt: {data['disk_attempt2']}",
            f"    # Disk reservation for 13th attempt: {data['disk_attempt13']}",
            f"    googlebatch_boot_disk_gb: {data['disk_function']}",
            f"    # Machine typde defined on best case: {data['googlebatch_machine_types'][0]}, {data['googlebatch_machine_prices'][0]} $ / hour",
            f"    # Machine typde defined on worst case: {data['googlebatch_machine_types'][-1]}, {data['googlebatch_machine_prices'][-1]} $ / hour",
            f"    googlebatch_machine_type: {data['googlebatch_machine_types'][1]}",
            f"    # Machine types are locked for google-batch-api",
            f"    googlebatch_image_family: 'batch-cos'",
            f"    # Disk types are locked for google-batch-api",
            f"    googlebatch_boot_disk_type: 'hyperdisk-balanced'",
            f"    # Image project are defined based on both keys above. See docs.",
            f"    googlebatch_image_project: 'batch-custom-image'",
        ]

    return text
    
def chose_google_machine_type(
    known: pandas.DataFrame,
    machine_types_path: str | pathlib.Path,
    googlebatch_boot_disk_type: str = "hyperdisk-balanced",
) -> str:
    """Chose the right machine type"""
    logger.info(f"Selecting Google machine types from {machine_types_path=}")
    # Thank you Nils Knieling !
    # https://gcloud-compute.com/europe-west9.html
    # Copy -> Paste -> to CSV
    # Loading machine types CSV
    machine_types = pandas.read_csv(machine_types_path)
    machine_types = machine_types[["Name", "vCPU", "Memory", "#", "Hour"]]

    # Filtering machines which cannot be used with batch or paris
    machine_types["vCPU"] = [int(str(i).split(" ")[0]) for i in machine_types.vCPU]
    machine_types = machine_types[machine_types["#"] == 3]
    allowed_machine_types_hyperdisk_balanced: tuple[str] = (
        "c3", "c3d", "c4", "c4a", "c4d", "n4", "n4a", "n4d", "h3", "m3", "g4", "a3", "a4"
    )
    allowed_machine_types_hyperdisk_extreme: tuple[str] = (
        "c3", "c3d", "c4", "n4", "m3"
    )
    allowed_machine_types_hyperdisk_balanced_nvidia: tuple[str] = (
        "n4", "n4a", "n4d"
    )
    # Keep only machines that accept hyperdisks.
    if googlebatch_boot_disk_type == "hyperdisk-balanced":
        machine_types = machine_types[machine_types["Name"].str.startswith(allowed_machine_types_hyperdisk_balanced)]
    elif googlebatch_boot_disk_type == "hyperdisk-extreme":
        machine_types = machine_types[machine_types["Name"].str.startswith(allowed_machine_types_hyperdisk_extreme)]

    # Sorting to always return cheapest mechine
    machine_types = machine_types.sort_values(by="Hour", ascending=True)

    # Actually assessing best machine type
    machine_lists: list[str | None] = [None, None, None]
    price_lists: list[str | None] = [None, None, None]
    # Best case
    cpu = known["threads"].min()
    mem_gb = known["memory_gb"].min()
    for machine in machine_types.itertuples():
        if (cpu <= machine.vCPU) and (mem_gb <= machine.Memory):
            machine_lists[0] = machine.Name
            price_lists[0] = machine.Hour
            break

    # Mean case
    cpu = known["threads"].mean()
    mem_gb = known["memory_gb"].mean() 
    for machine in machine_types.itertuples():
        if (cpu <= machine.vCPU) and (mem_gb <= machine.Memory):
            machine_lists[1] = machine.Name
            price_lists[1] = machine.Hour
            break

    # Worst case
    cpu = known["threads"].max()
    mem_gb = known["memory_gb"].max()
    for machine in machine_types.itertuples():
        if (cpu <= machine.vCPU) and (mem_gb <= machine.Memory):
            machine_lists[2] = machine.Name
            price_lists[2] = machine.Hour
            break

    return machine_lists, price_lists


def find_tsv(
    benchmark_dir: str | pathlib.Path,
) -> typing.Generator[pathlib.Path, None, None]:
    """Find all benchmark files"""
    if isinstance(benchmark_dir, str):
        benchmark_dir = pathlib.Path(benchmark_dir)
    logger.debug(f"Looking for benchmarks in {str(benchmark_dir.resolve())}")

    for content in benchmark_dir.iterdir():
        if content.is_dir():
            yield from find_tsv(content)
        else:
            with content.open("r") as content_stream:
                first_line = content_stream.readline()
                if first_line.startswith("s\th:m:s\t"):
                    yield content


def save_profile(resources: str, path: str | pathlib.Path) -> None:
    """save workflow-profile to a given path"""
    logger.debug(f"Saving profile to {path=}")
    if isinstance(path, str):
        path = pathlib.Path(path)

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as yaml_out_stream:
        yaml_out_stream.write("\n".join(resources))

def build_parser() -> argparse.ArgumentParser:
    """Build command line parser"""
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=rich_argparse.ArgumentDefaultsRichHelpFormatter,
    )
    parser.add_argument(
        "-b", "--benchmark",
        default=pathlib.Path("benchmarks"),
        type=pathlib.Path,
        help="Directory containing the benchmark data.",
    )
    parser.add_argument(
        "-m", "--machines",
        default=pathlib.Path("paris_machine_types.csv"),
        type=pathlib.Path,
        help="Path to CSV file listing available machines",
    )
    parser.add_argument(
        "-p", "--profile",
        type=pathlib.Path,
        default=pathlib.Path("workflow-profile-auto"),
        help="Output directory where profiles are written.",
    )
    parser.add_argument(
        "-d", "--google-disks",
        type=str,
        default="hyperdisk-balanced",
        help="The type of disks used with the google batch api",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="count",
        default=0,
        help="Increase verbosity (use -v for INFO, -vv for DEBUG).",
    )
    return parser


def setup_logging(verbosity: int) -> None:
    """set-up logging levels"""
    logger.remove()
    if verbosity >= 2:
        level = "DEBUG"
    elif verbosity == 1:
        level = "INFO"
    else:
        level = "WARNING"
    logger.add(sys.stderr, level=level, format="{time:DD-MM-YY:HH:mm:ss} | {level} | {message}")


def main(
    benchmark_dir: str | pathlib.Path = "benchmarks",
    machine_types_path: str | pathlib.Path = "paris_machine_types.csv",
    outdir: str | pathlib.Path = "workflow-profiles-auto",
    googlebatch_boot_disk_type: str = "hyperdisk-balanced",
) -> None:
    if isinstance(outdir, str):
        outdir = pathlib.Path(outdir)

    tsv_paths: list[pathlib.Path] = list(find_tsv(benchmark_dir=benchmark_dir))
    benchmarks_df: pandas.DataFrame = cat_dataframes(tsv_paths)
    benchmarks_df = process_extended_benchmark(benchmarks_df)
    outpath = outdir / "big_table.csv"
    outpath.parent.mkdir(parents=True, exist_ok=True)
    benchmarks_df.to_csv(outpath, index=False)

    for mode in ("local", "slurm", "google"):
        resources = ["set-resources:"]
        for rule in sorted(set(benchmarks_df.rule_name)):
            resources += rule_profile(
                known=benchmarks_df[benchmarks_df["rule_name"] == rule],
                mode=mode,
                machine_types_path=machine_types_path,
                googlebatch_boot_disk_type=googlebatch_boot_disk_type,
            )
        save_profile(resources, str(outdir / f"{mode}" / "config.yaml"))

if __name__ == "__main__":
    args = build_parser().parse_args()
    setup_logging(args.verbose)
    main(benchmark_dir=args.benchmark, machine_types_path=args.machines, outdir=args.profile, googlebatch_boot_disk_type=args.google_disks,)
    logger.success("Done.")
