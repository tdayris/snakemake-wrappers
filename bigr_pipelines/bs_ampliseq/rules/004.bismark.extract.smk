rule bismark_methylation_extractor:
    input:
        "bismark/align/{sample}.PE.bam",
    output:
        mbias_r1="data_output/Bismark/{sample}.M-bias_R1.png",
        mbias_r2="data_output/Bismark/{sample}.M-bias_R2.png",
        mbias_report=temp("bismark/meth/{sample}.M-bias.txt"),
        splitting_report=temp("bismark/meth/{sample}_PE_splitting_report.txt"),
        methylome_CpG_cov=protected("data_output/Bismark/{sample}.cpg.cov.gz"),
        methylome_CpG_mlevel_bedGraph=protected(
            "data_output/Bismark/{sample}.cpg.bedGraph.gz"
        ),
        read_base_meth_state_cpg="data_output/Bismark/CpG_context_{sample}.txt.gz",
        read_base_meth_state_chg="data_output/Bismark/CHG_context_{sample}.txt.gz",
        read_base_meth_state_chh="data_output/Bismark/CHH_context_{sample}.txt.gz",
    threads: 1
    resources:
        mem_mb=get_4gb_per_attempt,
        time_min=get_3h_per_attempt,
        tmpdir="tmp",
    log:
        "logs/bismark/meth/{sample}.log",
    params:
        output_dir="bismark/meth",
        extra=config["bismark"].get("extract", "--gzip --comprehensive --bedGraph"),
    wrapper:
        "bio/bismark/bismark_methylation_extractor"
