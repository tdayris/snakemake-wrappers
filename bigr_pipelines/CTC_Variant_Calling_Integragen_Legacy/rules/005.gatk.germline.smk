"""
gatk 3.7 haplotypecaller (on each normal sample) : java -Xmx8g -jar GenomeAnalysisTK.jar -T HaplotypeCaller -R hg38_basechr.fa -I {normal1}.cutadapt.sorted.rmmarkdup.bam -o {normal1}.cutadapt.sorted.rmmarkdup.hc.g.vcf.gz -ERC GVCF -bamout {normal1}.cutadapt.sorted.rmmarkdup.hc.bam
"""


rule gatk_haplotype_caller:
    input:
        bam="sambamba/markdup/{sample}.{status}.bam",
        fasta=config["ref"]["fasta"],
    output:
        vcf=temp("gatk/haplotypecaller/{sample}.{status}.g.vcf.gz"),
        bam="gatk/haplotypecaller/{sample}.baseline.bam",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 10,
        java_mem_gb=lambda wildcards, attempt: attempt * 9,
        time_min=lambda wildcards, attempt: attempt * 75,
        tmpdir="tmp",
    group:
        "baseline_wbc_calling"
    log:
        "logs/haplotypecaller/{sample}.{status}.log",
    params:
        "-ERC GVCF",
    conda:
        "envs/conda/gatk3.yaml"
    shell:
        "java -Xmx{resources.java_mem_gb}GB "
        "-jar GenomeAnalysisTK.jar "
        "-T HaplotypeCaller "
        "{params} "
        "-R {input.fasta} "
        "-I {input.bam} "
        "-O {output.vcf} "
        "-bamout {output.bam} "
        "> {log} 2>&1"


"""
gatk 3.7 genotypegvcfs : java -Xmx8g -jar GenomeAnalysisTK.jar -T GenotypeGVCFs -R hg38_basechr.fa -V {sample1_normal}.cutadapt.sorted.rmmarkdup.hc.g.vcf.gz -V {sample2_normal}.cutadapt.sorted.rmmarkdup.hc.g.vcf.gz -V {sample3_normal}.cutadapt.sorted.rmmarkdup.hc.g.vcf.gz -o {patient}_CTC_all.cutadapt.sorted.rmmarkdup.hc.vcf.gz
"""


rule gatk_genotype_gvcf:
    input:
        normal="gatk/haplotypecaller/{sample}.baseline.g.vcf.gz",
        wbc="gatk/haplotypecaller/{sample}.wbc.g.vcf.gz",
        fasta=config["ref"]["fasta"],
    output:
        temp("gatk/genotype_gvcf/baseline_wbc/{sample}.g.vcf.gz"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 10,
        java_mem_gb=lambda wildcards, attempt: attempt * 9,
        time_min=lambda wildcards, attempt: attempt * 75,
        tmpdir="tmp",
    group:
        "baseline_wbc_calling"
    log:
        "logs/gatk/genotype_gvcf/baseline_wbc/{sample}.log",
    params:
        "",
    conda:
        "envs/conda/gatk3.yaml"
    shell:
        "java -Xmx{resources.java_mem_gb}GB "
        "-jar GenomeAnalysisTK.jar "
        "-T GenotypeGVCFs "
        "{params} "
        "-R {input.fasta} "
        "-V {input.baseline} "
        "-V {input.wbc} "
        "-o {output} "
        "> {log} 2>&1"


"""
 gatk 3.7 selectvariants : java -Xmx8g -jar GenomeAnalysisTK.jar -T SelectVariants -R hg38_basechr.fa -V {patient}_CTC_all.cutadapt.sorted.rmmarkdup.hc.vcf.gz -selectType SNP -o {patient}_CTC_all.cutadapt.sorted.rmmarkdup.hc.snps.vcf
"""


rule gatk_select_variants:
    input:
        gvcf="gatk/genotype_gvcf/baseline_wbc/{sample}.g.vcf.gz",
        fasta=config["ref"]["fasta"],
    output:
        temp("gatk/select_variants/baseline_wbc/{sample}.g.vcf"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 10,
        java_mem_gb=lambda wildcards, attempt: attempt * 9,
        time_min=lambda wildcards, attempt: attempt * 75,
        tmpdir="tmp",
    group:
        "baseline_wbc_variant_filters"
    log:
        "logs/gatk/select_variants/baseline_wbc/{sample}.log",
    params:
        "-selectType SNP",
    conda:
        "envs/conda/gatk3.yaml"
    shell:
        "java -Xmx{resources.java_mem_gb}GB "
        "-jar GenomeAnalysisTK.jar "
        "-T SelectVariants "
        "{params} "
        "-R {input.fasta} "
        "-V {input.gvcf} "
        "-o {output} "
        "> {log} 2>&1"


"""
gatk 3.7 variantfiltration : java -Xmx8g -jar GenomeAnalysisTK.jar -T VariantFiltration -R hg38_basechr.fa -V {patient}_CTC_all.cutadapt.sorted.rmmarkdup.hc.snps.vcf --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filterName "my_snp_filter" -o {patient}_CTC_all.cutadapt.sorted.rmmarkdup.hc.snps.filtered.vcf
"""


rule gatk_variant_filtration:
    input:
        gvcf="gatk/select_variants/baseline_wbc/{sample}.g.vcf",
        fasta=config["ref"]["fasta"],
    output:
        temp("gatk/variant_filtration/baseline_wbc/{sample}.g.vcf"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 10,
        java_mem_gb=lambda wildcards, attempt: attempt * 9,
        time_min=lambda wildcards, attempt: attempt * 75,
        tmpdir="tmp",
    group:
        "baseline_wbc_variant_filters"
    log:
        "logs/gatk/variant_filtration/baseline_wbc/{sample}.log",
    params:
        "--filterExpression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0' "
        "--filterName 'custom_snp_filter'",
    conda:
        "envs/conda/gatk3.yaml"
    shell:
        "java -Xmx{resources.java_mem_gb}GB "
        "-jar GenomeAnalysisTK.jar "
        "-T VariantFiltration "
        "{params} "
        "-R {input.fasta} "
        "-V {input.gvcf} "
        "-o {output} "
        "> {log} 2>&1"
