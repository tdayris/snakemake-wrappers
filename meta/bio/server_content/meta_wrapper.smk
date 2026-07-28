rule curl_list_available_resources:
    """get the list of downloadable resources"""
    output:
        temp("<tmp>/curl_list_available_resources/{web_address}.xml"),
    log:
        "<log>/curl_list_available_resources/{web_address}.log",
    benchmark:
        "<benchmark>/curl_list_available_resources/{web_address}.tsv",
    threads: 1
    params:
        url="{web_address}",
        extra="--silent",
    conda:
        workflow.source_path("resources/curl.yaml"),
    shell:
        "curl {params.extra} {params.url:q} > {output:q} 2> {log:q}"


rule go_yq_format_to_csv:
    """easier parsing with snakemake"""
    input:
        "<tmp>/curl_list_available_resources/{web_address}.xml",
    output:
        temp("<tmp>/go_yq_format_to_csv/{web_address}.csv"),
    log:
        "<log>/go_yq_format_to_csv/{web_address}.log",
    benchmark:
        "<benchmark>/go_yq_format_to_csv/{web_address}.tsv",
    threads: 1
    params:
        subcommand="eval",
        extra="",
        expression=str(
            '.ListBucketResult.Contents[] | '
            '[ .Key, .LastModified, .Size, (.Owner.ID? // ""), '
            '(.ChecksumAlgorithm? // ""), (.ChecksumType? // "")] | '
            '@csv'
        ),
    wrapper:
        "v9.15.0/utils/go-yq"


rule add_header_to_csv:
    input:
        "<tmp>/go_yq_format_to_csv/{web_address}.csv",
    output:
        "<results>/available_data/{web_address}.csv",
    log:
        "<log>/available_data/{web_address}.log",
    benchmark:
        "<benchmark>/available_data/{web_address}.tsv",
    threads: 1
    params:
        subcommand="cat rows",
        extra="<( echo 'Path,Date,Size,Owner,ChecksumAlgorithm,ChecksumType' ) ",
    wrapper:
        "v3.4.0/utils/xsv"


rule bioconvert_csv_to_xls:
    input:
        "<results>/available_data/{web_address}.csv",
    output:
        "<results>/available_data/{web_address}.xls",
    log:
        "<log>/available_data/{web_address}.log",
    benchmark:
        "<benchmark>/available_data/{web_address}.tsv"
    threads: 1
    params:
        converter="csv2xls",
        extra="",
    wrapper:
        "v9.15.0/bio/bioconvert"
