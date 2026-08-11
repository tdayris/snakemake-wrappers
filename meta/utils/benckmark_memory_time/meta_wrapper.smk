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


def get_benchmarks_per_tool(
    wildcards: snakemake.iocontainers.Wildcards, 
    input_dir: str = input_benchmark,
) -> dict[str, list[str]]:
    """get the bechnmark files for a given tool"""
    return {
        "table": list(get_generator_of_benchmarks(f"{input_dir}/{wildcards.tool}")),
    }




def get_rule_report_per_tool(
    wildcards: snakemake.iocontainers.Wildcards,
    input_dir: str = input_benchmark, 
) -> dict[str, tuple[str]]:
    """get list of rule report per tool"""
    tool: str = str(wildcards.tool)
    return {
        tool: tuple(
            f"<temp>/go_yq_sum_input_file_size/{tool}.{rule_report}.csv"
            for rule_report in get_generator_of_benchmarks(
                benchmark_dir=f"{input_dir}/{tool}", 
                stem=True,
            )
        )
    }

tools_tpl: tuple[str] = tuple(str(child.stem) for child in pathlib.Path(input_benchmark).iterdir())
rule_reports_tpl: tuple[str] = tuple(get_generator_of_benchmarks(input_benchmark, stem=True))
print(tools_tpl)
print(rule_reports_tpl)

wildcard_constraints:
    tool=r"|".join(tools_tpl),
    rule_report=r"|".join(rule_reports_tpl),

rule xsv_format_benchmarks_to_csv:
    """format from tsv to csv"""
    input:
        table="<input_benchmark>/{tool}/{rule_report}.tsv",
    output:
        temp("<temp>/xsv_format_benchmarks_to_csv/{tool}.{rule_report}.csv"),
    log:
        "<logs>/xsv_format_benchmarks_to_csv/{tool}_{rule_report}.log",
    benchmark:
        "<benchmarks>/xsv/format_benchmarks_to_csv_{tool}_{rule_report}.tsv",
    threads: 1
    params:
        subcommand="fmt",
        extra="",
    wrapper:
        "v3.4.0/utils/xsv"


use rule xsv_format_benchmarks_to_csv as xsv_extract_size_time_memory with:
    """keep only columns of interest"""
    input:
        table="<temp>/xsv_format_benchmarks_to_csv/{tool}.{rule_report}.csv",
    output:
        temp("<temp>/xsv_extract_size_time_memory/{tool}.{rule_report}.csv"),
    log:
        "<logs>/xsv_extract_size_time_memory/{tool}.{rule_report}.log",
    benchmark:
        "<benchmarks>/xsv/extract_size_time_memory_{tool}.{rule_report}.tsv"
    params:
        extra="rule_name,s,'h:m:s',max_rss,input_size_mb",
        subcommand="select",


rule go_yq_sum_input_file_size:
    """parse json formatted input_file_size column"""
    input:
        "<temp>/xsv_extract_size_time_memory/{tool}.{rule_report}.csv"
    output:
        temp("<temp>/go_yq_sum_input_file_size/{tool}.{rule_report}.csv"),
    log:
        "<logs>/go_yq_sum_input_file_size/{tool}.{rule_report}.log",
    benchmark:
        "<benchmarks>/go_yq/sum_input_file_size_{tool}.{rule_report}.tsv"
    threads: 1
    params:
        extra="",
        subcommand="eval",
        expression=".[].input_size_mb = ( to_entries | .[].input_size_mb as $entry ireduce (0; . + ($entry.value | tonumber))) | .[].tool = \"{tool}\" ",
    wrapper:
        "v9.15.0/utils/go-yq"

use rule xsv_format_benchmarks_to_csv as xsv_cat_rules with:
    """aggregate requirements per tool"""
    input:
        unpack(get_rule_report_per_tool),
    output:
        "<results>/benchmark/{tool}/rules.csv",
    log:
        "<logs>/xsv_cat_rules/{tool}.log"
    benchmark:
        "<benchmarks>/xsv/cat_rules_{tool}.tsv"
    params:
        extra="",
        subcommand="cat rows", 


use rule xsv_format_benchmarks_to_csv as xsv_summarize_benchmarks with:
    """produce desciptive statistics over benchmarks"""
    input:
        "<results>/benchmark/{tool}/rules.csv",
    output:
        "<results>/benchmark/{tool}/summary.csv",
    log:
        "<logs>/xsv_summarize_benchmarks/{tool}.log",
    benchmark:
        "<benchmarks>/xsv/summarize_benchmarks_{tool}.tsv"
    params:
        extra="--select s,max_rss,input_size_mb",
        subcommand="stats",
        

use rule xsv_cat_rules as xsv_complete_benchmark with:
    """aggregate all tools into a single table"""
    input:
        table=expand(
            "<results>/benchmark/{tool}/rules.csv",
            tool=tools_tpl,
        ),
    output:
        "<results>/benchmark/benchmarks.csv",
    log:
        "<logs>/xsv_complete_benchmark.log",
    benchmark:
        "<benchmarks>/xsv/complete_benchmark.tsv"
 

rule xan_hist_time_per_tool:
    input:
        "<results>/benchmark/benchmarks.csv",
    output:
        temp("<temp>/xan_hist_time_per_tool.txt"),
    log:
        "<logs>/xan_hist_time_per_tool.log",
    benchmark:
        "<benchmarks>/xan_hist/time_per_tool.tsv",
    params:
        groupby="tool 'sum(s) as time'",
        hist="-v time -l tool --name 'Time used by tools' --unit ' seconds'",
    conda:
        "resources/xan.yaml"
    shell:
        "( xan groupby {params.groupby} {input:q} | "
        "  xan hist {params.hist} ) > {output:q} 2> {log:q}"


use rule xan_hist_time_per_tool as xan_hist_mb_per_tool with:
    log:
        "<logs>/xan_hist_time_per_tool.log",
    benchmark:
        "<benchmarks>/xan_hist/time_per_tool.tsv",
    params:
        groupby="tool 'sum(max_rss) as memory'",
        hist="-v memory -l tool --name 'RAM used by tools' --unit ' mb'",
    


