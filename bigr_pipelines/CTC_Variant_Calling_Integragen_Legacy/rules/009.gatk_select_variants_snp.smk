"""
 gatk 3.7 selectvariants : java -Xmx8g -jar GenomeAnalysisTK.jar -T SelectVariants -R hg38_basechr.fa -V {patient}_CTC_all.cutadapt.sorted.rmmarkdup.hc.vcf.gz -selectType SNP -o {patient}_CTC_all.cutadapt.sorted.rmmarkdup.hc.snps.vcf
"""


rule gatk_select_variants_snp_wbc_baseline:
    input:
        gvcf="gatk/genotype_gvcf/{sample}_{version}_{manip}.g.vcf.gz",
        fasta=config["ref"]["fasta"],
    output:
        temp("gatk/select_variants/{sample}_{version}_{manip}.g.vcf"),
    threads: 1
    resources:
        mem_mb=4 * 1024,
        java_mem_gb=4 * 1024,
        time_min=get_6h_per_attempt,
        tmpdir=tmp,
    log:
        "logs/gatk/select_variants/{sample}_{version}_{manip}.log",
    params:
        extra="-selectType SNP",
        jar="/mnt/beegfs/userdata/t_dayris/GATK3.7/devs/GATK/GenomeAnalysisTK.jar",
        tmp=tmp,
    conda:
        str(workflow_source_dir / "envs" / "gatk.yaml")
    shell:
        "gatk "
        "-Xmx{resources.java_mem_gb}M "
        "-Djava.io.tmpdir=\"{params.tmp}\" "
        "-T SelectVariants "
        "{params.extra} "
        "-R {input.fasta} "
        "-V {input.gvcf} "
        "-o {output} "
        "> {log} 2>&1"


rule gatk_select_variants_snp_ctc_baseline:
    input:
        gvcf="gatk/genotype_gvcf/{sample}_{version}_{manip}_{nb}.g.vcf.gz",
        fasta=config["ref"]["fasta"],
    output:
        temp("gatk/select_variants/{sample}_{version}_{manip}_{nb}.g.vcf"),
    threads: 1
    resources:
        mem_mb=4 * 1024,
        java_mem_gb=4 * 1024,
        time_min=get_6h_per_attempt,
        tmpdir=tmp,
    log:
        "logs/gatk/select_variants/{sample}_{version}_{manip}_{nb}.log",
    params:
        extra="-selectType SNP",
        jar="/mnt/beegfs/userdata/t_dayris/GATK3.7/devs/GATK/GenomeAnalysisTK.jar",
        tmp=tmp,
    conda:
        str(workflow_source_dir / "envs" / "gatk.yaml")
    shell:
        "gatk "
        "-Xmx{resources.java_mem_gb}M "
        "-Djava.io.tmpdir=\"{params.tmp}\" "
        "-T SelectVariants "
        "{params.extra} "
        "-R {input.fasta} "
        "-V {input.gvcf} "
        "-o {output} "
        "> {log} 2>&1"


rule gatk_select_variants_snp_baseline:
    input:
        gvcf="gatk/genotype_gvcf/{sample}.baseline.g.vcf.gz",
        fasta=config["ref"]["fasta"],
    output:
        temp("gatk/select_variants/{sample}.baseline.g.vcf"),
    threads: 1
    resources:
        mem_mb=4 * 1024,
        java_mem_gb=4 * 1024,
        time_min=get_6h_per_attempt,
        tmpdir=tmp,
    log:
        "logs/gatk/select_variants/{sample}.baseline.log",
    params:
        extra="-selectType SNP",
        jar="/mnt/beegfs/userdata/t_dayris/GATK3.7/devs/GATK/GenomeAnalysisTK.jar",
        tmp=tmp,
    conda:
        str(workflow_source_dir / "envs" / "gatk.yaml")
    shell:
        "gatk "
        "-Xmx{resources.java_mem_gb}M "
        "-Djava.io.tmpdir=\"{params.tmp}\" "
        "-T SelectVariants "
        "{params.extra} "
        "-R {input.fasta} "
        "-V {input.gvcf} "
        "-o {output} "
        "> {log} 2>&1"