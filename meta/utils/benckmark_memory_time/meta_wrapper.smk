import pathlib
import typing

input_benchmark: str = config["pathvars"]["input_benchmark"]

def get_generator_of_benchmarks(
    benchmark_dir: pathlib.Path | str,
    stem: bool = False,
) -> typing.Generator[str, None, None]:
    """recursively search for tsv (benchmark) files"""
    if isinstance(benchmark_dir, str):
        benchmark_dir = pathlib.Path(benchmark_dir)

    for child in benchmark_dir.iterdir():
        if child.is_dir():
            yield from get_generator_of_benchmarks(child, stem)
        elif str(child.name).endswith(".tsv"):
            yield str(child.stem) if stem else str(child)

data: list[str] = list(get_generator_of_benchmarks(input_benchmark))
print(data)


rule fd_search_tsv:
    input:
        input_benchmark,
    output:
        temp("<temp>/fd_search_tsv.txt"),
    log:
        "<logs>/fd_search_tsv.log",
    benchmark:
        "<benchmarks>/fd/search_tsv.tsv",
    threads: 7
    params:
        extra="",
    wrapper:
        "file:../../../../utils/fd"
        

rule xan_format_benchmarks_to_csv:
    input:
        data=data,
    output:
        temp("<temp>/xan_format_benchmarks_to_csv/benchmarks.csv"),
    log:
        "<logs>/xan_format_benchmarks_to_csv.log",
    benchmark:
        "<benchmarks>/xan/format_benchmarks_to_csv.tsv",
    threads: 7
    params:
        extra="--tee",
        expression=lambda wildcards, threads: str(
            f"cat rows --delimiter '\\t' --raw "
        ),
    wrapper:
        "file:../../../../utils/xan/run"
    

rule go_yq_sum_input_file_size:
    """parse json formatted input_file_size column"""
    input:
        "<temp>/xan_format_benchmarks_to_csv/benchmarks.csv",
    output:
        temp("<temp>/go_yq_sum_input_file_size/benchmarks.csv"),
    log:
        "<logs>/go_yq_sum_input_file_size/benchmarks.log",
    benchmark:
        "<benchmarks>/go_yq/sum_input_file_size_benchmarks.tsv"
    threads: 1
    params:
        extra="",
        subcommand="eval",
        expression='.[].input_size_mb = ( to_entries | .[].input_size_mb as $entry ireduce (0; . + ($entry.value | tonumber)))',
    wrapper:
        "v9.15.0/utils/go-yq"

use rule xan_format_benchmarks_to_csv as xan_compute_job_statistics with:
    """build max rss per io and max minutes per io"""
    input:
        data="<temp>/go_yq_sum_input_file_size/benchmarks.csv",
    output:
        "<results>/benchmarks/resources.csv",
    log:
        "<logs>/xan_compute_job_statistics.log",
    benchmark:
        "<benchmarks>/xan/compute_job_statistics.tsv",
    threads: 4
    params:
        extra="--tee",
        expression=str(
            "progress --title 'compute job statistics' | "
            "select --evaluate 'rule_name,max(0.5, s / 60) as minutes,max(0.5, input_size_mb) as input_size_mb,max(0.5, max_rss) as max_rss,s as seconds' | "
            "groupby rule_name 'max(minutes) as minutes,max(max_rss) as max_rss,max(input_size_mb) as input_size_mb' | "
            "select --evaluate 'rule_name,max_rss / input_size_mb as max_rss_per_mb,minutes / input_size_mb as minutes_per_mb,max_rss,minutes,input_size_mb' "
        ),


use rule xan_format_benchmarks_to_csv as xan_hist_rule with:
    """display statistics over max_rss and max_time"""
    input:
        data="<results>/benchmarks/resources.csv",
    output:
        temp("<temp>/xan_hist_rule/{rss_time}.txt"),
    log:
        "<logs>/xan_hist_rule/{rss_time}.log",
    benchmark:
        "<benchmarks>/xan_hist_rule/{rss_time}.tsv",
    wildcard_constraints:
        rss_time="max_rss|minutes",
    threads: 1
    params:
        extra="",
        expression="hist -v '{rss_time}' -l rule_name --name '{rss_time} per rule' --bar-size small --cols 120 --color 'always'",


use rule xan_format_benchmarks_to_csv as xan_display_table with:
    """print human readable table"""
    input:
        data="<results>/benchmarks/resources.csv",
    output:
        temp("<temp>/xan_display_table.csv"),
    log:
        "<logs>/xan_display_table.log",
    benchmark:
        "<benchmarks>/xan_display_table.tsv",
    threads: 1
    params:
        extra="",
        expression='view --all --name "Resource per MB of input" --select "rule_name,max_rss_per_mb,minutes_per_mb" --hide-index --color "always" --cols 120',

rule cat_xan_reports:
    input:
        expand(
            "<temp>/xan_hist_rule/{rss_time}.txt",
            rss_time=("max_rss", "minutes"),
        ),
        "<temp>/xan_display_table.csv",
    output:
        "<results>/benchmarks/resources.txt",
    log:
        "<logs>/cat_xan_reports.log",
    benchmark:
        "<benchmarks>/cat_xan_reports.tsv"
    threads: 1
    params:
        extra="",
    shell:
        "cat {params.extra} {input} > {output:q} 2> {log:q}"
    
        
