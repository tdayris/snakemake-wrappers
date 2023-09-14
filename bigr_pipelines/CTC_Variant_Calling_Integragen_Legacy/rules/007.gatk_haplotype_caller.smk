"""
gatk 3.7 haplotypecaller (on each normal sample) : java -Xmx8g -jar GenomeAnalysisTK.jar -T HaplotypeCaller -R hg38_basechr.fa -I {normal1}.cutadapt.sorted.rmmarkdup.bam -o {normal1}.cutadapt.sorted.rmmarkdup.hc.g.vcf.gz -ERC GVCF -bamout {normal1}.cutadapt.sorted.rmmarkdup.hc.bam
"""


rule gatk_haplotype_caller_ctc:
    input:
        bam="sambamba/markdup/{sample}_{version}_{manip}_{nb}.bam",
        fasta=config["ref"]["fasta"],
    output:
        vcf=temp("gatk/haplotypecaller/{sample}_{version}_{manip}_{nb}.g.vcf.gz"),
        bam="gatk/haplotypecaller/{sample}_{version}_{manip}_{nb}.bam",
    threads: 1
    resources:
        mem_mb=5 * 1024,
        java_mem_gb=4 * 1024,
        time_min=get_8h_per_attempt,
        tmpdir=tmp,
    log:
        "logs/haplotypecaller/{sample}_{version}_{manip}_{nb}.log",
    params:
        extra="-ERC GVCF",
        jar="/mnt/beegfs/userdata/t_dayris/GATK3.7/devs/GATK/GenomeAnalysisTK.jar",
        tmp=tmp,
    conda:
        str(workflow_source_dir / "envs" / "gatk.yaml")
    shell:
        "gatk "
        "-Xmx{resources.java_mem_gb}M "
        "-Djava.io.tmpdir=\"{params.tmp}\" "
        "-T HaplotypeCaller "
        "{params.extra} "
        "-R {input.fasta} "
        "-I {input.bam} "
        "-o {output.vcf} "
        "-bamout {output.bam} "
        "> {log} 2>&1"


rule gatk_haplotype_caller_wbc:
    input:
        bam="sambamba/markdup/{sample}_{version}_{manip}.bam",
        fasta=config["ref"]["fasta"],
    output:
        vcf=temp("gatk/haplotypecaller/{sample}_{version}_{manip}.g.vcf.gz"),
        bam="gatk/haplotypecaller/{sample}_{version}_{manip}.bam",
    threads: 1
    resources:
        mem_mb=5 * 1024,
        java_mem_gb=4 * 1024,
        time_min=get_8h_per_attempt,
        tmpdir=tmp,
    log:
        "logs/haplotypecaller/{sample}_{version}_{manip}.log",
    params:
        extra="-ERC GVCF",
        jar="/mnt/beegfs/userdata/t_dayris/GATK3.7/devs/GATK/GenomeAnalysisTK.jar",
        tmp=tmp,
    conda:
        str(workflow_source_dir / "envs" / "gatk.yaml")
    shell:
        "gatk "
        "-Xmx{resources.java_mem_gb}M "
        "-Djava.io.tmpdir=\"{params.tmp}\" "
        "-T HaplotypeCaller "
        "{params.extra} "
        "-R {input.fasta} "
        "-I {input.bam} "
        "-o {output.vcf} "
        "-bamout {output.bam} "
        "> {log} 2>&1"


rule gatk_haplotype_caller_baseline:
    input:
        bam="sambamba/markdup/{sample}.baseline.bam",
        fasta=config["ref"]["fasta"],
    output:
        vcf=temp("gatk/haplotypecaller/{sample}.baseline.g.vcf.gz"),
        bam="gatk/haplotypecaller/{sample}.baseline.bam",
    threads: 1
    resources:
        mem_mb=5 * 1024,
        java_mem_gb=4 * 1024,
        time_min=get_8h_per_attempt,
        tmpdir=tmp,
    log:
        "logs/haplotypecaller/{sample}.baseline.log",
    params:
        extra="-ERC GVCF",
        jar="/mnt/beegfs/userdata/t_dayris/GATK3.7/devs/GATK/GenomeAnalysisTK.jar",
        tmp=tmp,
    conda:
        str(workflow_source_dir / "envs" / "gatk.yaml")
    shell:
        "gatk "
        "-Xmx{resources.java_mem_gb}M "
        "-Djava.io.tmpdir=\"{params.tmp}\" "
        "-T HaplotypeCaller "
        "{params.extra} "
        "-R {input.fasta} "
        "-I {input.bam} "
        "-o {output.vcf} "
        "-bamout {output.bam} "
        "> {log} 2>&1"