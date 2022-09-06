rule bismark_reports:
    input:
        alignment_report=temp("bismark/bams/{sample}_report.txt"),
        nucleotide_report=temp("bismark/bams/{sample}.nucleotide_stats.txt"),
        mbias_report=temp("bismark/meth/{sample}.M-bias.txt"),
        splitting_report=temp("bismark/meth/{sample}_splitting_report.txt"),
    output:
        html="data_output/Bismark/{sample}.html",
    threads: 1
    resources:
        mem_mb=get_2gb_per_attempt,
        time_min=get_1h_per_attempt,
        tmpdir="tmp",
    params:
        skip_optional_reports=config["bismark"].get("skip_optional_reports", False),
    log:
        "logs/bismark/report/{sample}.log",
    wrapper:
        "bio/bismark/bismark2report"
