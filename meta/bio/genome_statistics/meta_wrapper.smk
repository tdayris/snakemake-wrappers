rule agat_sp_statistics:
    """gather genome-wide statistics"""
    input:
        gff="<genome_annotation>",
        gs="<genome_sequence>",
        config="<agat_config>",
    output:
        report="<reference>/<species>.<build>.<release>/statistics/statistics.txt",
        yaml=temp("<tmp>/agat_sp_statistics/<species>.<build>.<release>.yaml"),
        plot=directory("<reference>/<species>.<build>.<release>/statistics/graphs"),
    log:
        "<log>/agat_sp_statistics/<species>.<build>.<release>.log",
    benchmark:
        "<benchmark>/agat_sp_statistics/<species>.<build>.<release>.tsv",
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
        assembly_stats=temp("<tmp>/assembly_stats/<species>.<build>.<release>.tsv"),
    log:
        "<log>/assembly_stats/<species>.<build>.<release>.log",
    benchmark:
        "<benchmark>/assembly_stats/<species>.<build>.<release>.tsv",
    threads: 1
    params:
        extra="-t"
    wrapper:
        "v2.9.1/bio/assembly-stats"


rule go_yq_format_tsv_tp_yaml:
    """in order to include assembly stats in agat ones"""
    input:
        "<tmp>/assembly_stats/<species>.<build>.<release>.tsv",
    output:
        temp("<tmp>/go_yq_format_tsv/<species>.<build>.<release>.yaml"),
    log:
        "<log>/go_yq_format_tsv/<species>.<build>.<release>.log",
    benchmark:
        "<benchmark>/go_yq_format_tsv/<species>.<build>.<release>.tsv"
    threads: 1
    params:
        extra="",
        subcommand="",
        expression="",
    wrapper:
        "v9.14.0/utils/go-yq"


rule go_yq_include_assembly_into_agat:
    """in order to have a single entry genome statistics"""
    input:
        "<tmp>/agat_sp_statistics/<species>.<build>.<release>.yaml",
        "<tmp>/go_yq_format_tsv/<species>.<build>.<release>.yaml",
    output:
        "<reference>/<species>.<build>.<release>/statistics/statistics.yaml",
    default_target: True
    log:
        "<log>/go_yq_include_assembly_into_agat/<species>.<build>.<release>.log",
    benchmark:
        "<benchmark>/go_yq_include_assembly_into_agat/<species>.<build>.<release>.tsv"
    params:
        extra="--no-doc",
        subcommand="eval-all",
        expression='[.]',
    wrapper:
        "v9.14.0/utils/go-yq"
        
