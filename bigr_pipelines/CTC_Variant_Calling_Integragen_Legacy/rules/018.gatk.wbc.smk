rule gatk_select_variants_wbconly:
    input:
        vcf="grep/baseline_wbc/{sample}.vcf",
        fasta=config["ref"]["fasta"],
    output:
        temp("gatk/select_variants/wbc/{sample}.tmp.vcf"),
    threads: 1
    resources:
        mem_mb=get_2gb_per_attempt,
        time_min=get_45min_per_attempt,
        java_mem_gb=get_2gb_per_attempt,
        tmpdir=tmp,
    group:
        "retrieve_wbc"
    log:
        "logs/gatk/select_variants/wbc/{sample}.pass.log",
    params:
        extra="-sn {sample}.wbc",
        jar="/mnt/beegfs/userdata/t_dayris/GATK3.7/devs/GATK/GenomeAnalysisTK.jar",
        tmp=tmp,
    conda:
        str(workflow_source_dir / "envs" / "gatk.yaml")
    shell:
        # "java -Xmx{resources.java_mem_gb}GB "
        # "-jar {params.jar} "
        "gatk "
        "-Xmx{resources.java_mem_gb}GB "
        "-Djava.io.tmpdir=\"{params.tmp}\" "
        "-T SelectVariants "
        "{params.extra} "
        "-R {input.fasta} "
        "-V {input.vcf} "
        "-o {output} "
        "> {log} 2>&1"




rule grep_out_homozygote_wbc:
    input:
        "gatk/select_variants/wbc/{sample}.tmp.vcf",
    output:
        temp("gatk/select_variants/wbc/{sample}.vcf"),
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_35min_per_attempt,
        tmpdir=tmp,
    group:
        "retrieve_wbc"
    log:
        "logs/grep/{sample}.log",
    conda:
        str(workflow_source_dir / "envs" / "bash.yaml")
    params:
        ' -v "./.:0,0:0:0,0,0"',
    shell:
        "grep {params} {input} > {output} 2> {log}"


rule zip_wbc_variants:
    input:
        "gatk/select_variants/wbc/{sample}.vcf"
    output:
        protected("data_output/WBC/{sample}.vcf.gz")
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


rule tabix_wbc_variants:
    input:
        "data_output/wbc/{sample}.vcf.gz"
    output:
        protected("data_output/WBC/{sample}.vcf.gz.tbi")
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