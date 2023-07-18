"""
mutect 2.0 : java -Xmx8g -jar GenomeAnalysisTK.jar -T MuTect2 -I:tumor {sample_tumor}.cutadapt.sorted.bam -I:normal {sample_normal}.cutadapt.sorted.bam --output_mode EMIT_VARIANTS_ONLY -o {sample_tumor}vs{sample_normal}_mutect2.vcf.gz -R /mnt/beegfs/userdata/m_deloger/B20002_FRFA_01/DBs/bwa/hg38_basechr.fa --max_alt_alleles_in_normal_count 2 --max_alt_allele_in_normal_fraction 0.04 --maxReadsInRegionPerSample 100000
"""


rule mutect2:
    input:
        unpack(get_trio_wbc),
    output:
        "gatk/mutect2_wbc/{sample}.vcf.gz",
    threads: 1
    resources:
        mem_mb=4 * 1024,
        java_mem_gb=4 * 1024,
        time_min=get_6h_per_attempt,
        tmpdir=tmp,
    log:
        "logs/gatk/mutect2_wbc/{sample}.log",
    params:
        extra=(
            "--max_alt_alleles_in_normal_count 2 "
            "--max_alt_allele_in_normal_fraction 0.04 "
            "--maxReadsInRegionPerSample 100000 "
            "--output_mode EMIT_VARIANTS_ONLY "
        ),
        tmp=tmp,
    conda:
        str(workflow_source_dir / "envs" / "gatk.yaml")
    shell:
        "gatk "
        "-Xmx{resources.java_mem_gb}M "
        "-Djava.io.tmpdir=\"{params.tmp}\" "
        "{params.extra} "
        "-T MuTect2 "
        "-I:tumor {input.tumor} "
        "-I:normal {input.normal} "
        "-o {output} "
        "-R {input.fasta} "
        "> {log} 2>&1 "


rule tabix_mutect2:
    input:
        "gatk/mutect2_wbc/{sample}.vcf.gz",
    output:
        "gatk/mutect2_wbc/{sample}.vcf.gz.tbi",
    threads: 1
    resources:
        mem_mb=get_4gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir=tmp,
    log:
        "logs/tabix/{sample}.mutect2_wbc.log",
    params:
        "-p vcf",
    wrapper:
        "bio/tabix/index"


rule unzip_mutect2:
    input:
        "gatk/mutect2_wbc/{sample}.vcf.gz",
    output:
        temp("gatk/mutect2_wbc/{sample}.vcf"),
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_35min_per_attempt,
        tmpdir=tmp,
    conda:
        str(workflow_source_dir / "envs" / "bash.yaml")
    log:
        "logs/gatk/mutect2_wbc/unzip/{sample}.log",
    params:
        "--decompress --force --verbose --stdout",
    shell:
        "gzip {params} {input} > {output} 2> {log}"
