"""
ensembl VEP 87.0 refseq (on each normal VCF/TSV) : singularity run -B /mnt/beegfs:/mnt/beegfs /mnt/beegfs/software/vep/87/ensembl-vep_release_87.0.sif vep -i {patient}_CTC_all.cutadapt.sorted.rmmarkdup.hc.snps.filtered.PASS.{sample1_normal}.vcf -o {patient}_CTC_all.cutadapt.sorted.rmmarkdup.hc.snps.filtered.PASS.{sample1_normal}.annotated.vcf --species homo_sapiens --assembly GRCh38 --cache --dir_cache /mnt/beegfs/database/bioinfo/vep/87/ --everything --offline --fasta /mnt/beegfs/database/bioinfo/vep/87/homo_sapiens/87_GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa --vcf --refseq
"""


rule ensembl_vep_bam_readcount_ctc:
    input:
        vcf="brc2vep/filtered/{sample}_{version}_{manip}_{nb}.tsv",
        cache=config["ref"]["vep"],
        fasta="resources/GRCh38.fasta",
    output:
        vcf=temp("vep/annotate/{sample}_{version}_{manip}_{nb}.brc.vcf"),
    threads: 1
    resources:
        mem_mb=get_10gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir=tmp,
    log:
        "logs/vep/{sample}_{version}_{manip}_{nb}.log",
    params:
        "--species homo_sapiens "
        "--assembly GRCh38 "
        "--cache "
        "--everything "
        "--offline "
        "--vcf "
        "--force_overwrite "
        "--format ensembl "
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


rule ensembl_vep_bam_readcount_wbc:
    input:
        vcf="brc2vep/filtered/{sample}_{version}_{manip}.tsv",
        cache=config["ref"]["vep"],
        fasta="resources/GRCh38.fasta",
    output:
        vcf=temp("vep/annotate/{sample}_{version}_{manip}.brc.vcf"),
    threads: 1
    resources:
        mem_mb=get_10gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir=tmp,
    log:
        "logs/vep/{sample}_{version}_{manip}.log",
    params:
        "--species homo_sapiens "
        "--assembly GRCh38 "
        "--cache "
        "--everything "
        "--offline "
        "--vcf "
        "--force_overwrite "
        "--format ensembl "
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


rule ensembl_vep_bam_readcount_baseline:
    input:
        vcf="brc2vep/filtered/{sample}.baseline.tsv",
        cache=config["ref"]["vep"],
        fasta="resources/GRCh38.fasta",
    output:
        vcf=temp("vep/annotate/{sample}.baseline.brc.vcf"),
    threads: 1
    resources:
        mem_mb=get_10gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir=tmp,
    log:
        "logs/vep/{sample}.baseline.log",
    params:
        "--species homo_sapiens "
        "--assembly GRCh38 "
        "--cache "
        "--everything "
        "--offline "
        "--vcf "
        "--force_overwrite "
        "--format ensembl "
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