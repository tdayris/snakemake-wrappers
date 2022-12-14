rule salmon_quantification:
    input:
        **unpack(get_salmon_index),
        **unpack(get_salmon_fastq),
    output:
        quant="030.quantification/salmon_{genome_build}.{genome_release}/{sample}/quant.sf",
        quant_genes="030.quantification/salmon_{genome_build}.{genome_release}/{sample}/quant.genes.sf",
        lib_format="030.quantification/salmon_{genome_build}.{genome_release}/{sample}/lib_format_counts.json",
        cmd_info="030.quantification/salmon_{genome_build}.{genome_release}/{sample}/cmd_info.json",
        logs="030.quantification/salmon_{genome_build}.{genome_release}/{sample}/logs/salmon_quant.log",
        flendist="030.quantification/salmon_{genome_build}.{genome_release}/{sample}/libParams/flenDist.txt",
        ambig_info="030.quantification/salmon_{genome_build}.{genome_release}/{sample}/aux_info/ambig_info.tsv",
        bootstraps=directory("030.quantification/salmon_{genome_build}.{genome_release}/{sample}/aux_info/bootstrap"),
        exp=multiext("030.quantification/salmon_{genome_build}.{genome_release}/{sample}/aux_info/exp", "3_pos.gz", "3_seq.gz", "5_pos.gz", "5_seq.gz", "cted_bias.gz", "_gc.gz"),
        obs=multiext("030.quantification/salmon_{genome_build}.{genome_release}/{sample}/aux_info/", "3_pos.gz", "3_seq.gz", "5_pos.gz", "5_seq.gz", "erved_bias_3p.gz", "erved_bias.gz", "_gc.gz"),
        expected_bias="030.quantification/salmon_{genome_build}.{genome_release}/{sample}/aux_info/expected_bias.gz",
        fld="030.quantification/salmon_{genome_build}.{genome_release}/{sample}/aux_info/fld.gz",
        meta_info="030.quantification/salmon_{genome_build}.{genome_release}/{sample}/aux_info/meta_info.json",
    threads: 20
    resources:
        mem_mb=get_10gb_and_10gb_per_attempt,
        time_min=get_45min_and_20min_per_attempt,
        tmpdir="tmp"
    log:
        "logs/030.quantification.salmon/{sample}.{genome_build}.{genome_release}.smk"
    params:
        extra="--numBootstraps 100 --gcBias --seqBias --posBias"
    wrapper:
        "v1.20.0/bio/salmon/quant"



            salmon_logs=snakemake.input["salmon_logs"],
            salmon_metas=snakemake.input["salmon_meta"],
            salmon_mapping_plot=snakemake.output["mapping_status_mqc"],
            salmon_assigned_fragments=snakemake.output["mapping_counts_mqc"],
            salmon_mapping_rates_general_stat=snakemake.output["salmon_general_table"],


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

    