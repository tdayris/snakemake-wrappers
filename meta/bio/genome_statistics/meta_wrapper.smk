rule agat_sp_statistics:
    """gather genome-wide statistics"""
    input:
        gff="<genome_annotation>",
        gs="<genome_sequence>",
        config="<agat_config>",
    output:
        report="<resources>/<species>.<build>.<release>/statistics/statistics.txt",
        yaml=temp("<temp>/agat_sp_statistics/<species>.<build>.<release>.yaml"),
        plot=directory("<resources>/<species>.<build>.<release>/statistics/graphs"),
    log:
        "<logs>/agat_sp_statistics/<species>.<build>.<release>.log",
    benchmark:
        "<benchmarks>/agat/sp_statistics_<species>.<build>.<release>.tsv"
    threads: 6
    params:
        command="agat_sp_statistics.pl",
    wrapper:
        "v9.6.0/bio/agat"


rule assembly_stats:
    """gather contigs-wide statistics"""
    input:
        assembly="<genome_sequence>",
    output:
        assembly_stats=temp("<temp>/assembly_stats/<species>.<build>.<release>.tsv"),
    log:
        "<logs>/assembly_stats/<species>.<build>.<release>.log",
    benchmark:
        "<benchmarks>/assembly_stats/<species>.<build>.<release>.tsv"
    threads: 1
    params:
        extra="-t",
    wrapper:
        "v2.9.1/bio/assembly-stats"


rule go_yq_format_tsv_to_yaml:
    """in order to include assembly stats in agat ones"""
    input:
        "<temp>/assembly_stats/<species>.<build>.<release>.tsv",
    output:
        temp("<temp>/go_yq_format_tsv/<species>.<build>.<release>.yaml"),
    log:
        "<logs>/go_yq_format_tsv/<species>.<build>.<release>.log",
    benchmark:
        "<benchmarks>/go_yq/format_tsv_<species>.<build>.<release>.tsv"
    threads: 1
    params:
        extra="",
        subcommand="",
        expression="",
    wrapper:
        "v9.14.0/utils/go-yq"


use rule go_yq_format_tsv_to_yaml as go_yq_include_assembly_into_agat with:
    """in order to have a single entry genome statistics"""
    input:
        "<temp>/agat_sp_statistics/<species>.<build>.<release>.yaml",
        "<temp>/go_yq_format_tsv/<species>.<build>.<release>.yaml",
    output:
        "<resources>/<species>.<build>.<release>/statistics/statistics.yaml",
    log:
        "<logs>/go_yq_include_assembly_into_agat/<species>.<build>.<release>.log",
    benchmark:
        "<benchmarks>/go_yq/include_assembly_into_agat_<species>.<build>.<release>.tsv"
    params:
        extra="--no-doc",
        subcommand="eval-all",
        expression="[.]",
