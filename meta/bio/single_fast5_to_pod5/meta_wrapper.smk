import csv
import pandas

# If no samples are provided in name space, then
# load samples table with column separator auto-detection
if not (("samples_table" in locals()) or ("samples_table" in globals())):
    samples_table_path: str = config.get("samples_table", "config/samples.csv")

    # Detect separator
    with open(samples_table_path, "r") as table_stream:
        dialect: csv.Dialect = csv.Sniffer().sniff(table_stream.readline())
        table_stream.seek(0)

    # Load samples table
    samples_table: pandas.DataFrame = pandas.read_csv(
        samples_table_path,
        sep=dialect.delimiter,
        header=0,
        index_col=None,
        comment="#",
        dtype=str,
    )

# If there are any, list **all** required columns missing from samples table
required_tables: set[str] = {"sample_id", "fast5_archive"}
missing_cols: set[str] = required_tables - set(samples_table.columns)
if len(missing_cols) > 0:
    raise KeyError(f"Could not find {missing_cols=} in {samples_table.columns=}")
samples_table.set_index("sample_id", inplace=True, drop=False)


rule single_fast5_to_pod5_aria_download_single_read_archive:
    """let fast5 data be available for later use"""
    output:
        temp("temp/access_fast5_archive/{sample}.tar.gz"),
    log:
        "logs/single_fast5_to_pod5_aria_download_single_read_archive/{sample}.log",
    benchmark:
        "benchmarks/aria2/single_fast5_to_pod5_aria_download_single_read_archive_{sample}.tsv"
    params:
        url=lambda wildcards: samples_table.at[wildcards.sample, "fast5_archive"],
        extra="--retry-wait 5",
    wrapper:
        f"{config['repo']}/utils/aria2c"


rule single_fast5_to_pod5_untar_single_read_archive:
    """decompress all single read fast5 samples"""
    input:
        archive="temp/access_fast5_archive/{sample}.tar.gz",
    output:
        outdir=temp(directory("temp/single_fast5_to_pod5_untar_single_read_archive/{sample}")),
    log:
        "log/single_fast5_to_pod5_untar_single_read_archive/{sample}.log",
    benchmark:
        "benchmarks/tar/single_fast5_to_pod5_untar_single_read_archive_{sample}.tsv"
    params:
        extra="--extract --overwrite --verbose --force-local",
        explicit_owner=True,
        explicit_group=True,
    wrapper:
        f"{config['repo']}/utils/tar"


rule single_fast5_to_pod5_fast5_aggregation_with_ont_api:
    """reformat single read fast5 into multi read fast5"""
    input:
        indir="temp/single_fast5_to_pod5_untar_single_read_archive/{sample}",
    output:
        outdir=temp(directory("temp/single_fast5_to_pod5_fast5_aggregation_with_ont_api/{sample}")),
    log:
        "logs/single_fast5_to_pod5_fast5_aggregation_with_ont_api/{sample}.log",
    benchmark:
        "benchmarks/single_to_multi_fast5/single_fast5_to_pod5_fast5_aggregation_with_ont_api/{sample}.tsv"
    container:
        # "docker://quay.io/biocontainers/ont-fast5-api:0.4.1--py35_0"
        "docker://nanozoo/ont-fast5-api:3.1.6--a980386"
    threads: 15
    params:
        extra="--recursive --batch_size 1000000",
    shell:
        "single_to_multi_fast5 {params.extra} "
        "--input_path {input.indir:q} "
        "--save_path {output.outdir:q} "
        "--threads {threads} "
        "> {log:q} 2>&1 "


rule single_fast5_to_pod5_fast5_to_pod5_conversion:
    """reformat mutl-read fast5 into pod5"""
    input:
        indir="temp/single_fast5_to_pod5_fast5_aggregation_with_ont_api/{sample}",
    output:
        pod=temp("temp/single_fast5_to_pod5_fast5_to_pod5_conversion/{sample}.pod5"),
    log:
        "logs/single_fast5_to_pod5_fast5_to_pod5_conversion/{sample}.log",
    benchmark:
        "benchmarks/pod5/single_fast5_to_pod5_fast5_to_pod5_conversion/{sample}.tsv"
    container:
        #"docker://quay.io/biocontainers/pod5:0.3.44--pyhdfd78af_0"
        "docker://staphb/pod5:0.3.49"
    threads: 15
    params:
        extra="--recursive --force-overwrite ",
        mk="--parents --verbose",
        outdir=lambda wildcards, output: subpath(output.pod, parent=True),
    shell:
        "mkdir {params.mk} {params.outdir:q} > {log:q} 2>&1 && "
        "pod5 convert from_fast5 {params.extra} --threads {threads} "
        "--output {output.pod:q} {input.indir:q} >> {log:q} 2>&1 "


rule single_fast5_to_pod5_material_and_methods_single_fast5_to_pod5:
    """provide information on single_fast5_to_pod5 conversion methods"""
    input:
        rst="resources/reports/material_and_methods_single_fast5_to_pod5.rst",
    output:
        html="resources/{species}.{build}.{release}/sequences/Single_fast5_to_pod5_material_and_methods.html",
    log:
        "logs/single_fast5_to_pod5_material_and_methods_single_fast5_to_pod5/{species}.{build}.{release}.log",
    benchmark:
        "benchmarks/rst2html/single_fast5_to_pod5_material_and_methods_single_fast5_to_pod5_{species}.{build}.{release}.tsv"
    group:
        "mapping_dorado_material_and_methods"
    threads: 1
    params:
        extra="--title='Material and methods for dorado mapping' --embed-stylesheet",
    wrapper:
        f"{config['repo']}/utils/docutils/rst2html"
