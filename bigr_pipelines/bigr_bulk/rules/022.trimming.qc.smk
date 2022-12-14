rule multiqc_trimming:
    input:
        expand("020.trimming/fastp/html/{sample}.{ext}", sample=sample_list, ext=["html", "json"])
    output:
        "data_output/MultiQC/Trimming.html",
        directory("data_output/MultiQC/Trimming_data"),
    threads: 1
    resources:
        mem_mb=get_1gb_per_gb_of_input,
        time_min=get_10min_per_go,
        tmpdir="tmp"
    log:
        "logs/022.trimming.qc.log"
    params:
        extra="-f",
        use_input_files_only=True,
    wrapper:
        "v1.20.0/bio/multiqc"