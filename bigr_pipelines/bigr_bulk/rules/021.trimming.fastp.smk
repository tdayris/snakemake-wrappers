rule fastp_trimming_single_ended:
    input:
        sample=["data_input/{sample}.fq.gz"]
    output:
        trimmed="020.trimming/fastp/se/{sample}.fastq",
        html="020.trimming/fastp/html/{sample}.html",
        json="020.trimming/fastp/json/{sample}.json",
    threads: 20
    resources:
        mem_mb=get_4go_per_attempt,
        time_min=get_10min_per_go,
        tmpdir="tmp",
    log:
        "logs/021.fastp/{sample}.log"
    params:
        adapters=None,
        extra=(
            "--cut_front "
            "--cut_tail "
            "--cut_window_size 6 "
            "--cut_mean_quality 10 "
            "--unqualified_percent_limit 50 "
            "--n_base_limit 7 "
            "--average_qual 0 "
            "--length_required 15 "
            "--overrepresentation_analysis"
        ),
    wrapper:
        "v1.20.0/bio/fastp"


use rule fastp_trimming_single_ended as fastp_trimming_pair_ended with:
    input:
        sample=expand(
            "data_input/{sample}.{stream}.fq.gz",
            stream=stream_list,
            allow_missing=True,
        ),
    output:
        trimmed=expand(
            "020.trimming/fastp/pe/{sample}.{stream}.fastq",
            stream=stream_list,
            allow_missing=True,
        ),
        html="020.trimming/fastp/html/{sample}.html",
        json="020.trimming/fastp/json/{sample}.json",