
rule ensembl_vep_haplotype_caller_baseline:
    input:
        vcf="data_output/Baseline/{sample}.vcf.gz",
        vcf_tbi="data_output/Baseline/{sample}.vcf.gz.tbi",
        cache=config["ref"]["vep"],
        fasta="resources/GRCh38.fasta",
    output:
        vcf=temp("vep/annotate/{sample}.baseline.vcf"),
    threads: 1
    resources:
        mem_mb=get_20gb_per_attempt,
        time_min=get_3h_per_attempt,
        tmpdir=tmp,
    log:
        "logs/vep/{sample}.log",
    params:
        "--species homo_sapiens "
        "--assembly GRCh38 "
        "--cache "
        "--everything "
        "--offline "
        "--vcf "
        "--format vcf "
        "--force_overwrite "
        "--refseq ",
    container:
        "/mnt/beegfs/software/vep/87/ensembl-vep_release_87.0.sif"
    shell:
        "vep {params} "
        "--input_file {input.vcf} "
        "--output_file {output.vcf} "
        "--dir_cache {input.cache} "
        "--fasta {input.fasta} "
        "> {log} 2>&1 "


rule ensembl_vep_haplotype_caller_wbc:
    input:
        vcf="data_output/WBC/{sample}.vcf.gz",
        vcf_tbi="data_output/WBC/{sample}.vcf.gz.tbi",
        cache=config["ref"]["vep"],
        fasta="resources/GRCh38.fasta",
    output:
        vcf=temp("vep/annotate/{sample}.wbc.vcf"),
    threads: 1
    resources:
        mem_mb=get_20gb_per_attempt,
        time_min=get_3h_per_attempt,
        tmpdir=tmp,
    log:
        "logs/vep/{sample}.log",
    params:
        "--species homo_sapiens "
        "--assembly GRCh38 "
        "--cache "
        "--everything "
        "--offline "
        "--vcf "
        "--format vcf "
        "--force_overwrite "
        "--refseq ",
    container:
        "/mnt/beegfs/software/vep/87/ensembl-vep_release_87.0.sif"
    shell:
        "vep {params} "
        "--input_file {input.vcf} "
        "--output_file {output.vcf} "
        "--dir_cache {input.cache} "
        "--fasta {input.fasta} "
        "> {log} 2>&1 "