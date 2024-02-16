"""
gatk 3.7 variantfiltration : java -Xmx8g -jar GenomeAnalysisTK.jar -T VariantFiltration -R hg38_basechr.fa -V {patient}_CTC_all.cutadapt.sorted.rmmarkdup.hc.snps.vcf --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filterName "my_snp_filter" -o {patient}_CTC_all.cutadapt.sorted.rmmarkdup.hc.snps.filtered.vcf
"""


rule gatk_variant_filtration_wbc_baseline:
    input:
        gvcf="gatk/genotype_gvcf/{sample}_{version}_{manip}.g.vcf.gz",
        fasta=config["ref"]["fasta"],
    output:
        temp("gatk/variant_filtration/{sample}_{version}_{manip}.g.vcf"),
    threads: 1
    resources:
        mem_mb=4 * 1024,
        java_mem_gb=4 * 1024,
        time_min=get_6h_per_attempt,
        tmpdir=tmp,
    log:
        "logs/gatk/variant_filtration/{sample}_{version}_{manip}.log",
    params:
        extra=(
            "--filterExpression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0' "
            "--filterName 'custom_snp_filter'"
        ),
        jar="/mnt/beegfs/userdata/t_dayris/GATK3.7/devs/GATK/GenomeAnalysisTK.jar",
        tmp=tmp,
    conda:
        str(workflow_source_dir / "envs" / "gatk.yaml")
    shell:
        "gatk "
        "-Xmx{resources.java_mem_gb}M "
        "-Djava.io.tmpdir=\"{params.tmp}\" "
        "-T VariantFiltration "
        "{params.extra} "
        "-R {input.fasta} "
        "-V {input.gvcf} "
        "-o {output} "
        "> {log} 2>&1"



rule gatk_variant_filtration_ctc_baseline:
    input:
        gvcf="gatk/genotype_gvcf/{sample}_{version}_{manip}_{nb}.g.vcf.gz",
        fasta=config["ref"]["fasta"],
    output:
        temp("gatk/variant_filtration/{sample}_{version}_{manip}_{nb}.g.vcf"),
    threads: 1
    resources:
        mem_mb=4 * 1024,
        java_mem_gb=4 * 1024,
        time_min=get_6h_per_attempt,
        tmpdir=tmp,
    log:
        "logs/gatk/variant_filtration/{sample}_{version}_{manip}_{nb}.log",
    params:
        extra=(
            "--filterExpression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0' "
            "--filterName 'custom_snp_filter'"
        ),
        jar="/mnt/beegfs/userdata/t_dayris/GATK3.7/devs/GATK/GenomeAnalysisTK.jar",
        tmp=tmp,
    conda:
        str(workflow_source_dir / "envs" / "gatk.yaml")
    shell:
        "gatk "
        "-Xmx{resources.java_mem_gb}M "
        "-Djava.io.tmpdir=\"{params.tmp}\" "
        "-T VariantFiltration "
        "{params.extra} "
        "-R {input.fasta} "
        "-V {input.gvcf} "
        "-o {output} "
        "> {log} 2>&1"


rule gatk_variant_filtration_baseline:
    input:
        gvcf="gatk/genotype_gvcf/{sample}.baseline.g.vcf.gz",
        fasta=config["ref"]["fasta"],
    output:
        temp("gatk/variant_filtration/{sample}.baseline.g.vcf"),
    threads: 1
    resources:
        mem_mb=4 * 1024,
        java_mem_gb=4 * 1024,
        time_min=get_6h_per_attempt,
        tmpdir=tmp,
    log:
        "logs/gatk/variant_filtration/{sample}.baseline.log",
    params:
        extra=(
            "--filterExpression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0' "
            "--filterName 'custom_snp_filter'"
        ),
        jar="/mnt/beegfs/userdata/t_dayris/GATK3.7/devs/GATK/GenomeAnalysisTK.jar",
        tmp=tmp,
    conda:
        str(workflow_source_dir / "envs" / "gatk.yaml")
    shell:
        "gatk "
        "-Xmx{resources.java_mem_gb}M "
        "-Djava.io.tmpdir=\"{params.tmp}\" "
        "-T VariantFiltration "
        "{params.extra} "
        "-R {input.fasta} "
        "-V {input.gvcf} "
        "-o {output} "
        "> {log} 2>&1"
