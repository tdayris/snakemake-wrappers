rule quality_control_results:
    input:
        fastp=expand(
            "002.fastp/{ext}/{sample}.fastp.{ext}",
            sample=sample_list,
            ext=fastp_ext,
        ),
        fastq_screen=expand(
            "003.fastq_screen/{sample}.{stream}.fastq_screen.{ext}",
            sample=sample_list,
            stream=streams,
            ext=fastqscreen_ext,
        ),
    output:
        "data_output/multiqc/MultiQC.QC.html",
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir="tmp",
    log:
        "logs/snakefile.multiqc/salmon.log",
    wrapper:
        "bio/multiqc"
