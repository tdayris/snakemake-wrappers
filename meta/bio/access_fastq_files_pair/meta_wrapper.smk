import collections.abc
import csv
import os.path
import pandas
import snakemake.iocontainers
from snakemake_interface_storage_plugins.storage_provider import StorageProviderBase

# Define user-readable types
PandasIndexValueType = collections.abc.Hashable
PandasColumnValueType = collections.abc.Hashable
StorageProviderType = StorageProviderBase


# Load samples table with column separator auto-detection
samples_table_path: str = config.get("samples_table", "config/samples.csv")
with open(samples_table_path, "r") as table_stream:
    dialect: csv.Dialect = csv.Sniffer().sniff(table_stream.readline())
    table_stream.seek(0)

samples_table: pandas.DataFrame = pandas.read_csv(
    samples_table_path,
    header=0,
    index_col=None,
    sep=dialect.delimiter,
    comment="#",
    dtype=str,
)

# If there are any, list **all** required columns missing from samples table
required_tables: set[str] = {"sample_id", "upstream_file", "downstream_file"}
missing_cols: set[str] = required_tables - set(samples_table.columns)
if len(missing_cols) > 0:
    raise KeyError(f"Could not find {missing_cols=} in {samples_table.columns=}")
samples_table.set_index("sample_id", inplace=True)
    

# Building possible storage access methods with plugins
default_irods_env: str = os.path.expanduser("~/.irods/irods_environment.json")

#storage irods_files:
#    provider="irods",
#    storage_irods_host=config.get("irods_host", default_irods_env),

storage local_files:
    provider="fs",

storage http_files:
    provider="http",

def get_provider(path: str) -> StorageProviderType:
    """This function select the correct storage provider"""
    if path.startswith("http"):
        return storage.http_files(path)
    if path.startswith(tuple(config["irods_prefix"])):
        return storage.irods_files(path)
    if path.startswith(tuple(config["cold_storage"])):
        return storage.local_files(path)
    return path


# Functions to access samples fastq files within samples table
def get_fastq_from_pandas(
    index: PandasIndexValueType,
    column: PandasColumnValueType,
    samples_table: pandas.DataFrame = samples_table,
) -> StorageProviderType | list[StorageProviderType]:
    """This function parses the sample table"""
    fastq = samples_table.at[index, column]
    if "," in fastq:
        return [get_provider(p) for p in fastq.split(",")]
    return [get_provider(fastq)]


def get_fastq_rule_inteface(
    wildcards: snakemake.iocontainers.Wildcards,
    samples_table: pandas.DataFrame = samples_table,
) -> StorageProviderType | list[StorageProviderType]:
    """
    This function is an interface for rule unpacking.
    It lists the input files with their storage providers.
    """
    stream = str(wildcards.stream)
    return get_fastq_from_pandas(
        index = str(wildcards.sample),
        column = str(
            "upstream_file" if int(wildcards.stream) == 1 else "downstream_file"
        ),
        samples_table = samples_table,
    )


# Actual rules
rule fastcat_fastq:
    """concatenate resequencings into a single file"""
    input:
        unpack(get_fastq_rule_inteface),
    output:
        fastq=temp("<tmp>/fastcat_fastq/{sample}.{stream}.fastq.gz"),
        basecallers=temp(
            "<tmp>/fastcat_fastq/fastcat_{sample}.{stream}/basecallers.tsv"
        ),
        per_file=temp(
            "<tmp>/fastcat_fastq/fastcat_{sample}.{stream}/per_file.tsv"
        ),
        runids=temp(
            "<tmp>/fastcat_fastq/fastcat_{sample}.{stream}/runids.tsv"
        ),
        quality=temp(
            "<tmp>/fastcat_fastq/fastcat_{sample}.{stream}/histograms/quality.hist"
        ),
        length=temp(
            "<tmp>/fastcat_fastq/fastcat_{sample}.{stream}/histograms/length.hist"
        ),
    log:
        "<log>/fastcat_fastq/{sample}.{stream}.log",
    benchmark:
        "<benchmark>/fastcat_fastq/{sample}.{stream}.tsv"
    threads: 2
    conda:
        workflow.source_path("resources/fastcat.yaml"),
    params:
        extra=lambda wildcards: str(
            "--input-fmt 'fastq' --dust "
            f"--sample '{wildcards.sample}.{wildcards.stream}' "
        ),
        crabz_extra="--format 'gzip'",
        temp_out=lambda wildcards, output: str(Path(str(output.basecallers)).parent),
    shell:
        # Compression starts when fastcat process is over
        # so there is no threads over-use
        "( rm -rf {params.temp_out:q} && "
        "  fastcat fastq {params.extra} --threads {threads} "
        "  --output {params.temp_out:q} {input} | "
        "  crabz {params.crabz_extra} --output {output[0]:q} "
        "  --compression-threads {threads} & "
        " ) > {log} 2>&1"



rule xsv_cat_fastcat_per_file:
    """aggregate fastcat per sample qc table"""
    input:
        table=expand(
            "<tmp>/fastcat_fastq/fastcat_{sample}.{stream}/per_file.tsv",
            sample=samples_table.index,
            stream={1, 2},
        ),
    output:
        "<results>/fastq_qc/fastcat_per_fastq_file_qc.csv",
    log:
        "<log>/xsv_cat_fastcat_per_file.log",
    benchmark:
        "<benchmark>/xsv_cat_fastcat_per_file.tsv"
    threads: 1
    params:
        subcommand="cat rows",
        extra="",
    wrapper:
        "v3.4.0/utils/xsv"


rule xsv_frequency_fastcat_per_file:
    """summarize the fastcat aggregated table"""
    input:
        table="<results>/fastq_qc/fastcat_per_fastq_file_qc.csv",
    output:
        "<results>/fastq_qc/fastcat_per_fastq_file_qc_stats.csv",
    log:
        "<log>/xsv_frequency_fastcat_per_file.log",
    benchmark:
        "<benchmark>/xsv_frequency_fastcat_per_file.tsv",
    threads: 1
    params:
        subcommand="stats",
        extra="--select 'n_seqs,n_bases,min_length,max_length,mean_quality,f_file_ok,f_stream_error,f_qual_missing,f_qual_truncated,f_unknown_error,r_record_seen,r_record_ok,r_too_long,r_too_short,r_low_quality,r_dust_masked'",
    wrapper:
        "v3.4.0/utils/xsv"
