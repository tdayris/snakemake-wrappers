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


rule aria_download_single_read_archive:
    """let fast5 data be available for later use"""
    output:
       temp("<temp>/access_fast5_archive/<sample>.tar.gz"),
    log:
        "<logs>/access_fast5_archive/<sample>.log",
    benchmark:
        "<benchmarks>/aria2/access_fast5_archive_<sample>.tsv"
    params:
        url=lambda wildcards: samples_table.at[wildcards.sample, "fast5_archive"],
        extra="--retry-wait 5",
    wrapper:
        "v6.2.0/utils/aria2c"


rule untar_single_read_archive:
    """decompress all single-read fast5 samples"""
    input:
        "<temp>/access_fast5_archive/<sample>.tar.gz",
    output:
        temp(directory("<temp>/untar_single_read_archive/<sample>")),
    log:
        "log/untar_single_read_archive/<sample>.log"
    benchmark:
        "<benchmarks>/untar_single_read_archive/<sample>.tsv"
    params:
        extra=str(
            "--no-auto-compress --extract --overwrite --gzip "
            "--verbose --no-same-owner --no-same-permissions --force-local"
        ),
    shell:
        "tar {params.extra} --file {input:q} --directory {output:q} > {log:q} 2>&1"


rule fast5_aggregation_with_ont_api:
    """reformat single-read fast5 into multi-read fast5"""
    input:
        indir="<temp>/untar_single_read_archive/<sample>",
    output:
        outdir=temp(directory("<temp>/fast5_aggregation_with_ont_api/<sample>")),
    log:
        "<logs>/fast5_aggregation_with_ont_api/<sample>.log",
    benchmark:
        "<benchmarks>/fast5_aggregation_with_ont_api/<sample>.tsv"
    container:
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


rule fast5_to_pod5_conversion:
    """reformat mutl-read fast5 into pod5"""
    input:
        indir="<temp>/fast5_aggregation_with_ont_api/<sample>",
    output:
        pod=temp("<temp>/fast5_to_pod5_conversion/<sample>.pod5"),
    log:
        "<logs>/fast5_to_pod5_conversion/<sample>.log",
    benchmark:
        "<benchmarks>/fast5_to_pod5_conversion/<sample>.tsv"
    container:
        "docker://staphb/pod5:0.3.36"
    threads: 15
    params:
        extra="--recursive --force-overwrite ",
        mk="--parents --verbose",
        outdir=lambda wildcards, output: subpath(output.pod, parent=True),
    shell:
        "mkdir {params.mk} {params.outdir:q} > {log:q} 2>&1 && "
        "pod5 convert from_fast5 {params.extra} --threads {threads} "
        "--output {output.pod:q} {input.indir:q} >> {log:q} 2>&1 "
