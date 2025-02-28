"""
gatk 3.7 genotypegvcfs : java -Xmx8g -jar GenomeAnalysisTK.jar -T GenotypeGVCFs -R hg38_basechr.fa -V {sample1_normal}.cutadapt.sorted.rmmarkdup.hc.g.vcf.gz -V {sample2_normal}.cutadapt.sorted.rmmarkdup.hc.g.vcf.gz -V {sample3_normal}.cutadapt.sorted.rmmarkdup.hc.g.vcf.gz -o {patient}_CTC_all.cutadapt.sorted.rmmarkdup.hc.vcf.gz
"""


rule gatk_genotype_gvcf_baseline:
    input:
        normal="gatk/haplotypecaller/{sample}.baseline.g.vcf.gz",
        fasta=config["ref"]["fasta"],
    output:
        temp("gatk/genotype_gvcf/{sample}.baseline.g.vcf.gz"),
    threads: 1
    resources:
        mem_mb=4 * 1024,
        java_mem_gb=4 * 1024,
        time_min=get_6h_per_attempt,
        tmpdir=tmp,
    log:
        "logs/gatk/genotype_gvcf/{sample}.baseline.log",
    params:
        extra="",
        jar="/mnt/beegfs/userdata/t_dayris/GATK3.7/devs/GATK/GenomeAnalysisTK.jar",
        tmp=tmp,
    conda:
        str(workflow_source_dir / "envs" / "gatk.yaml")
    shell:
        "gatk "
        "-Xmx{resources.java_mem_gb}M "
        "-Djava.io.tmpdir=\"{params.tmp}\" "
        "-T GenotypeGVCFs "
        "{params.extra} "
        "-R {input.fasta} "
        "-V {input.normal} "
        "-o {output} "
        "> {log} 2>&1"


rule gatk_genotype_gvcf_wbc_baseline:
    input:
        normal="gatk/haplotypecaller/{sample}.baseline.g.vcf.gz",
        wbc="gatk/haplotypecaller/{sample}_{version}_{manip}.g.vcf.gz",
        fasta=config["ref"]["fasta"],
    output:
        temp("gatk/genotype_gvcf/{sample}_{version}_{manip}.g.vcf.gz"),
    threads: 1
    resources:
        mem_mb=4 * 1024,
        java_mem_gb=4 * 1024,
        time_min=get_6h_per_attempt,
        tmpdir=tmp,
    log:
        "logs/gatk/genotype_gvcf/{sample}_{version}_{manip}.log",
    params:
        extra="",
        jar="/mnt/beegfs/userdata/t_dayris/GATK3.7/devs/GATK/GenomeAnalysisTK.jar",
        tmp=tmp,
    conda:
        str(workflow_source_dir / "envs" / "gatk.yaml")
    shell:
        "gatk "
        "-Xmx{resources.java_mem_gb}M "
        "-Djava.io.tmpdir=\"{params.tmp}\" "
        "-T GenotypeGVCFs "
        "{params.extra} "
        "-R {input.fasta} "
        "-V {input.normal} "
        "-V {input.wbc} "
        "-o {output} "
        "> {log} 2>&1"


rule gatk_genotype_gvcf_ctc_baseline:
    input:
        normal="gatk/haplotypecaller/{sample}.baseline.g.vcf.gz",
        ctc="gatk/haplotypecaller/{sample}_{version}_{manip}_{nb}.g.vcf.gz",
        fasta=config["ref"]["fasta"],
    output:
        temp("gatk/genotype_gvcf/{sample}_{version}_{manip}_{nb}.g.vcf.gz"),
    threads: 1
    resources:
        mem_mb=4 * 1024,
        java_mem_gb=4 * 1024,
        time_min=get_6h_per_attempt,
        tmpdir=tmp,
    log:
        "logs/gatk/genotype_gvcf/{sample}_{version}_{manip}_{nb}.log",
    params:
        extra="",
        jar="/mnt/beegfs/userdata/t_dayris/GATK3.7/devs/GATK/GenomeAnalysisTK.jar",
        tmp=tmp,
    conda:
        str(workflow_source_dir / "envs" / "gatk.yaml")
    shell:
        "gatk "
        "-Xmx{resources.java_mem_gb}M "
        "-Djava.io.tmpdir=\"{params.tmp}\" "
        "-T GenotypeGVCFs "
        "{params.extra} "
        "-R {input.fasta} "
        "-V {input.normal} "
        "-V {input.ctc} "
        "-o {output} "
        "> {log} 2>&1"