rule salmon_to_multiqc:
    input:
        salmon_logs=expand("030.quantification/salmon_{genome_build}.{genome_release}/{sample}/lib_format_counts.json", sample=rna_sample_list, allow_missing=True),
        salmon_meta=expand("030.quantification/salmon_{genome_build}.{genome_release}/{sample}/aux_info/meta_info.json", sample=rna_sample_list, allow_missing=True),
    output:
        mapping_status_mqc=temp("multiqc_config/{genome_build}.{genome_release}/mapping_status_mqc.json"),
        mapping_counts_mqc=temp("multiqc_config/{genome_build}.{genome_release}/mapping_counts_mqc.json"),
        salmon_general_table=temp("multiqc_config/{genome_build}.{genome_release}/general_tabl_mqc.json"),
    threads: 1
    resources:
        mem_mb=get_3gb_per_gb,
        time_min=get_10min_per_go,
        tmpdir="tmp",
    log:
        "logs/030.quantification.salmon/{sample}.{genome_build}.{genome_release}_mqc_config.log"
    conda:
        str(workflow_source_dir / "envs" / "030.python.yaml")
    script:
        str(workflow_source_dir / "scripts" / "030.multiqc.salmon.log.py")


rule multiqc_salmon:
    input:
        expand("020.trimming/fastp/html/{sample}.{ext}", sample=sample_list, ext=["html", "json"]),
        expand("multiqc_config/{genome_build}.{genome_release}/mapping_status_mqc.json"),
        expand("multiqc_config/{genome_build}.{genome_release}/mapping_counts_mqc.json"),
        expand("multiqc_config/{genome_build}.{genome_release}/general_tabl_mqc.json"),
    output:
        protected("data_output/MultiQC/Quantification.{genome_build}.{genome_release}.html"),
        protected(directory("data_output/MultiQC/Quantification.{genome_build}.{genome_release}_data")),
    threads: 1
    resources:
        mem_mb=get_1gb_per_gb_of_input,
        time_min=get_10min_per_go,
        tmpdir="tmp"
    log:
        "logs/022.trimming.qc.log"
    params:
        extra="--force --module fastp --module salmon --title 'RNA-Seq quantification quality control report'",
        use_input_files_only=True,
    wrapper:
        "v1.20.0/bio/multiqc"