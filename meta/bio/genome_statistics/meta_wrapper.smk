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
        f"{config['repo']}/bio/agat"


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
        f"{config['repo']}/bio/assembly-stats"


rule go_yq_format_tsv_to_yaml:
    """in order to include assembly stats in agat ones"""
    input:
        "<temp>/{rule_content}/<species>.<build>.<release>.tsv",
    output:
        temp("<temp>/go_yq_format_tsv/<species>.<build>.<release>.{rule_content}.yaml"),
    log:
        "<logs>/go_yq_format_tsv/<species>.<build>.<release>.{rule_content}.log",
    benchmark:
        "<benchmarks>/go_yq/format_tsv_<species>.<build>.<release>.{rule_content}.tsv"
    threads: 1
    params:
        extra="",
        subcommand="",
        expression="",
    wrapper:
        f"{config['repo']}/utils/go-yq"


rule seqkit_stats_on_genome_sequence:
    """get detailed information on the genome sequence with seqkit"""
    input:
        fastx="<genome_sequence>",
    output:
        stats=temp("<temp>/seqkit_stats_on_genome_sequence/<species>.<build>.<release>.tsv"),
    log:
        "<logs>/seqkit_stats_on_genome_sequence/<species>.<build>.<release>.log",
    benchmark:
        "<benchmarks>/seqkit_stats/genome_sequence.<species>.<build>.<release>.tsv"
    threads: 5
    params:
        command="stats",
        extra="--all --tabular",
    wrapper:
        f"{config['repo']}/bio/seqkit"


rule fastaguard_on_genome_sequence:
    """get detailed information on the genome sequence with fastaguard"""
    input:
        fasta="<genome_sequence>",
    output:
        html="<resources>/<species>.<build>.<release>/statistics/fastaguard_report.html",
        json="<resources>/<species>.<build>.<release>/statistics/fastaguard_report.json",
        tsv=temp("<temp>/fastaguard_on_genome_sequence/<species>.<build>.<release>.tsv"),
        multiqc=temp("<temp>/fastaguard_on_genome_sequence/<species>.<build>.<release>._mqc.json"),
    log:
        "<logs>/fastaguard_on_genome_sequence/<species>.<build>.<release>.log",
    benchmark:
        "<benchmarks>/fastaguard/genome_sequence.<species>.<build>.<release>.tsv",
    params:
        extra="--profile assembly --gate pipeline",
    wrapper:
        f"{config['repo']}/bio/fastaguard"


use rule go_yq_format_tsv_to_yaml as go_yq_include_assembly_into_agat with:
    """in order to have a single entry genome statistics"""
    input:
        "<temp>/agat_sp_statistics/<species>.<build>.<release>.yaml",
        "<temp>/go_yq_format_tsv/<species>.<build>.<release>.assembly_stats.yaml",
        "<temp>/go_yq_format_tsv/<species>.<build>.<release>.seqkit_stats_on_genome_sequence.yaml",
        "<temp>/go_yq_format_tsv/<species>.<build>.<release>.fastaguard_on_genome_sequence.yaml",
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
