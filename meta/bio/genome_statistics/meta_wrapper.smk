def local_lookup(dpath, within, default):
    """Return local(path) if path is within config"""
    tmp = lookup(dpath=dpath, within=within, default=None)
    if tmp is None:
        return default
    return local(tmp)


rule genome_statistics_agat_sp_statistics:
    """gather genome-wide statistics"""
    input:
        gff=local_lookup(
            dpath="pathvars/genome_annotation",
            within=config,
            default="resources/{species}.{build}.{release}/annotations/{species}.{build}.{release}.gff3",
        ),
        gs=local_lookup(
            dpath="pathvars/genome_sequence",
            within=config,
            default="resources/{species}.{build}.{release}/sequences/{species}.{build}.{release}.dna.fasta",
        ),
        config=local_lookup(
            dpath="pathvars/agat_config",
            within=config,
            default="resources/{species}.{build}.{release}/annotations/{species}.{build}.{release}.agat_gff3_config.yaml",
        ),
    output:
        report="resources/{species}.{build}.{release}/statistics/statistics.txt",
        yaml=temp("temp/genome_statistics_agat_sp_statistics/{species}.{build}.{release}.yaml"),
        plot=directory("resources/{species}.{build}.{release}/statistics/graphs"),
    log:
        "logs/genome_statistics_agat_sp_statistics/{species}.{build}.{release}.log",
    benchmark:
        "benchmarks/agat/genome_statistics_agat_sp_statistics_{species}.{build}.{release}.tsv"
    threads: 6
    params:
        command="agat_sp_statistics.pl",
    wrapper:
        f"{config['repo']}/bio/agat"


rule genome_statistics_assembly_stats:
    """gather contigs-wide statistics"""
    input:
        assembly=local_lookup(
            dpath="pathvars/genome_sequence",
            within=config,
            default="resources/{species}.{build}.{release}/sequences/{species}.{build}.{release}.dna.fasta",
        ),
    output:
        assembly_stats=temp("temp/genome_statistics_assembly_stats/{species}.{build}.{release}.tsv"),
    log:
        "logs/genome_statistics_assembly_stats/{species}.{build}.{release}.log",
    benchmark:
        "benchmarks/assembly_stats/genome_statistics_assembly_stats/{species}.{build}.{release}.tsv"
    threads: 1
    params:
        extra="-t",
    wrapper:
        f"{config['repo']}/bio/assembly-stats"


rule genome_statistics_go_yq_format_tsv_to_yaml:
    """in order to include assembly stats in agat ones"""
    input:
        "temp/{rule_content}/{species}.{build}.{release}.tsv",
    output:
        temp("temp/go_yq_format_tsv/{species}.{build}.{release}.{rule_content}.yaml"),
    log:
        "logs/go_yq_format_tsv/{species}.{build}.{release}.{rule_content}.log",
    benchmark:
        "benchmarks/go_yq/genome_statistics_go_yq_format_tsv_to_yaml_{species}.{build}.{release}.{rule_content}.tsv"
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
        fastx=local_lookup(
            dpath="pathvars/genome_sequence",
            within=config,
            default="resources/{species}.{build}.{release}/sequences/{species}.{build}.{release}.dna.fasta",
        ),
    output:
        stats=temp(
            "temp/seqkit_stats_on_genome_sequence/{species}.{build}.{release}.tsv"
        ),
    log:
        "logs/seqkit_stats_on_genome_sequence/{species}.{build}.{release}.log",
    benchmark:
        "benchmarks/seqkit_stats/seqkit_stats_on_genome_sequence.{species}.{build}.{release}.tsv"
    threads: 5
    params:
        command="stats",
        extra="--all --tabular",
    wrapper:
        f"{config['repo']}/bio/seqkit"


rule genome_statistics_fastaguard_on_genome_sequence:
    """get detailed information on the genome sequence with fastaguard"""
    input:
        fasta=local_lookup(
            dpath="pathvars/genome_sequence",
            within=config,
            default="resources/{species}.{build}.{release}/sequences/{species}.{build}.{release}.dna.fasta",
        ),
    output:
        html=report(
            "resources/{species}.{build}.{release}/statistics/fastaguard_report.html",
            caption=workflow.source_path("resources/reports/fastaguard_report.rst"),
            category="Genome information",
            subcategory="{species}.{build}.{release}",
        ),
        json="resources/{species}.{build}.{release}/statistics/fastaguard_report.json",
        tsv=temp("temp/genome_statistics_fastaguard_on_genome_sequence/{species}.{build}.{release}.tsv"),
        multiqc=temp(
            "temp/genome_statistics_fastaguard_on_genome_sequence/{species}.{build}.{release}._mqc.json"
        ),
    log:
        "logs/genome_statistics_fastaguard_on_genome_sequence/{species}.{build}.{release}.log",
    benchmark:
        "benchmarks/fastaguard/genome_statistics_fastaguard_on_genome_sequence.{species}.{build}.{release}.tsv"
    group:
        "fastaguard_qc"
    threads: 1
    params:
        extra="--profile assembly --gate pipeline",
    wrapper:
        f"{config['repo']}/bio/fastaguard"


use rule genome_statistics_go_yq_format_tsv_to_yaml as genome_statistics_go_yq_include_assembly_into_agat with:
    """in order to have a single entry genome statistics"""
    input:
        "temp/genome_statistics_agat_sp_statistics/{species}.{build}.{release}.yaml",
        "temp/go_yq_format_tsv/{species}.{build}.{release}.genome_statistics_assembly_stats.yaml",
        "temp/go_yq_format_tsv/{species}.{build}.{release}.seqkit_stats_on_genome_sequence.yaml",
        "temp/go_yq_format_tsv/{species}.{build}.{release}.genome_statistics_fastaguard_on_genome_sequence.yaml",
    output:
        report(
            "resources/{species}.{build}.{release}/statistics/statistics.yaml",
            caption=workflow.source_path("resources/reports/genome_statistics.rst"),
            category="Genome information",
            subcategory="{species}.{build}.{release}",
        ),
    log:
        "logs/genome_statistics_go_yq_include_assembly_into_agat/{species}.{build}.{release}.log",
    benchmark:
        "benchmarks/go_yq/genome_statistics_go_yq_include_assembly_into_agat_{species}.{build}.{release}.tsv"
    params:
        extra="--no-doc",
        subcommand="eval-all",
        expression="[.]",


rule genome_statistics_multiqc_genome_qualities:
    """report genome qualities from fastaguard with multiqc"""
    input:
        "resources/{species}.{build}.{release}/statistics/fastaguard_report.html",
    output:
        report(
            "resources/{species}.{build}.{release}/statistics/Multiqc.html",
            caption=workflow.source_path(
                "resources/reports/multiqc_genome_statistics.rst"
            ),
            category="Genome information",
            subcategory="{species}.{build}.{release}",
        ),
        "resources/{species}.{build}.{release}/statistics/Multiqc_data.zip",
    log:
        "logs/genome_statistics_multiqc_genome_qualities/{species}.{build}.{release}.log",
    benchmark:
        "benchmarks/multiqc/genome_statistics_multiqc_genome_qualities_{species}.{build}.{release}.tsv"
    group:
        "fastaguard_qc"
    threads: 1
    params:
        extra="--verbose",
        use_input_files_only=True,
    wrapper:
        f"{config['repo']}/bio/multiqc"
