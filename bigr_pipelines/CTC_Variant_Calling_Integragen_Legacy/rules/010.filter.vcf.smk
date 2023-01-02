rule gatk_genotype_gvcf_ctc:
    input:
        ctc="gatk/haplotypecaller/{sample}.ctc.g.vcf.gz",
        fasta=config["ref"]["fasta"],
    output:
        temp("gatk/genotype_gvcf/hc_ctc/{sample}.g.vcf.gz"),
    threads: 1
    resources:
        mem_mb=get_10gb_per_attempt,
        java_mem_gb=get_10gb_per_attempt,
        time_min=get_2h_per_attempt,
        tmpdir="tmp",
    log:
        "logs/gatk/genotype_gvcf/baseline_ctc/{sample}.log",
    params:
        extra="",
        jar="/mnt/beegfs/userdata/t_dayris/GATK3.7/devs/GATK/GenomeAnalysisTK.jar"
    conda:
        str(workflow_source_dir / "envs" / "gatk.yaml")
    shell:
        # "java -Xmx{resources.java_mem_gb}MB "
        # "-jar {params.jar} "
        "gatk "
        "-Xmx{resources.java_mem_gb}MB "
        "-Djava.io.tmpdir='{resources.tmpdir}' "
        "-T GenotypeGVCFs "
        "{params.extra} "
        "-R {input.fasta} "
        "-V {input.ctc} "
        "-o {output} "
        "> {log} 2>&1"


rule filter_haplotype_ctc_vcf_non_snp:
    input:
        gvcf="gatk/genotype_gvcf/hc_ctc/{sample}.g.vcf.gz",
        fasta=config["ref"]["fasta"],
    output:
        temp("gatk/select_variants/hc_ctc/{sample}.g.vcf"),
    threads: 1
    resources:
        mem_mb=get_10gb_per_attempt,
        java_mem_gb=get_10gb_per_attempt,
        time_min=get_2h_per_attempt,
        tmpdir="tmp",
    log:
        "logs/gatk/select_variants/baseline_ctc/{sample}.log",
    params:
        extra="-selectType SNP",
        jar="/mnt/beegfs/userdata/t_dayris/GATK3.7/devs/GATK/GenomeAnalysisTK.jar"
    conda:
        str(workflow_source_dir / "envs" / "gatk.yaml")
    shell:
        # "java -Xmx{resources.java_mem_gb}MB "
        # "-jar {params.jar} "
        "gatk "
        "-Xmx{resources.java_mem_gb}MB "
        "-Djava.io.tmpdir='{resources.tmpdir}' "
        "-T SelectVariants "
        "{params.extra} "
        "-R {input.fasta} "
        "-V {input.gvcf} "
        "-o {output} "
        "> {log} 2>&1"


rule filter_haplotype_ctc_vcf_custom:
    input:
        gvcf="gatk/select_variants/hc_ctc/{sample}.g.vcf",
        fasta=config["ref"]["fasta"],
    output:
        vcf=temp("gatk/select_variants/hc_ctc/{sample}.tmp.vcf"),
    threads: 1
    resources:
        mem_mb=get_10gb_per_attempt,
        java_mem_gb=get_10gb_per_attempt,
        time_min=get_2h_per_attempt,
        tmpdir="tmp",
    log:
        "logs/gatk/variant_filtration/hc_ctc/{sample}.log",
    params:
        extra=(
            "--filterExpression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0' "
            "--filterName 'custom_snp_filter'"
        ),
    conda:
        str(workflow_source_dir / "envs" / "gatk.yaml")
    shell:
        "gatk "
        "-Xmx{resources.java_mem_gb}MB "
        "-Djava.io.tmpdir='{resources.tmpdir}' "
        "-T VariantFiltration "
        "{params.extra} "
        "-R {input.fasta} "
        "-V {input.gvcf} "
        "-o {output} "
        "> {log} 2>&1"


rule grep_out_homozygote_ctc:
    input:
        "gatk/select_variants/hc_ctc/{sample}.tmp.vcf",
    output:
        pipe("gatk/select_variants/hc_ctc/{sample}.no_homoz.vcf"),
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_35min_per_attempt,
        tmpdir="tmp",
    log:
        "logs/grep/{sample}.log",
    conda:
        str(workflow_source_dir / "envs" / "bash.yaml")
    params:
        ' -v "./.:0,0:0:0,0,0"',
    shell:
        "grep {params} {input} > {output} 2> {log}"


rule grep_out_filtered_ctc:
    input:
        "gatk/select_variants/hc_ctc/{sample}.no_homoz.vcf"
    output:
        temp("gatk/select_variants/hc_ctc/{sample}.vcf")
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_35min_per_attempt,
        tmpdir="tmp",
    log:
        "logs/grep/{sample}.log",
    conda:
        str(workflow_source_dir / "envs" / "bash.yaml")
    params:
        '{if ($0 ~ "^#") {print $0} else {if ($7 == ".") {print $0} else {if ($7 == "PASS") {print $0} } } }',
    shell:
        "awk {params} {input} > {output} 2> {log}"


rule zip_variants_ctc:
    input:
        "gatk/select_variants/hc_ctc/{sample}.vcf"
    output:
        protected("data_output/HC_CTC/{sample}.vcf.gz")
    threads: 2
    resources:
        mem_mb=get_4gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir="tmp"
    log:
        "logs/bcftools/view/{sample}.baseline.log"
    params:
        extra=""
    wrapper:
        "bio/bcftools/view"


rule tabix_variants_ctc:
    input:
        "data_output/HC_CTC/{sample}.vcf.gz"
    output:
        protected("data_output/HC_CTC/{sample}.vcf.gz.tbi")
    threads: 1
    resources:
        mem_mb=get_4gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir="tmp"
    log:
        "logs/tabix/{sample}.baseline.log"
    params:
        "-p vcf"
    wrapper:
        "bio/tabix/index"