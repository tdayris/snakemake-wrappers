"""
mutect 2.0 : java -Xmx8g -jar GenomeAnalysisTK.jar -T MuTect2 -I:tumor {sample_tumor}.cutadapt.sorted.bam -I:normal {sample_normal}.cutadapt.sorted.bam --output_mode EMIT_VARIANTS_ONLY -o {sample_tumor}vs{sample_normal}_mutect2.vcf.gz -R /mnt/beegfs/userdata/m_deloger/B20002_FRFA_01/DBs/bwa/hg38_basechr.fa --max_alt_alleles_in_normal_count 2 --max_alt_allele_in_normal_fraction 0.04 --maxReadsInRegionPerSample 100000
"""


rule gatk_mutect2_ctc:
    input:
        tumor="sambamba/markdup/{sample}_{version}_{manip}_{nb}.bam",
        normal="sambamba/markdup/{sample}.baseline.bam",
        fastaconfig["ref"]["fasta"],
    output:
        "gatk/mutect2/{sample}_{version}_{manip}_{nb}.vcf.gz",
    threads: 1
    resources:
        mem_mb=4 * 1024,
        java_mem_gb=4 * 1024,
        time_min=get_6h_per_attempt,
        tmpdir=tmp,
    log:
        "logs/gatk/mutect2/{sample}_{version}_{manip}_{nb}.log",
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



rule gatk_mutect2_wbc:
    input:
        tumor="sambamba/markdup/{sample}_{version}_{manip}.bam",
        normal="sambamba/markdup/{sample}.baseline.bam",
        fastaconfig["ref"]["fasta"],
    output:
        "gatk/mutect2/{sample}_{version}_{manip}.vcf.gz",
    threads: 1
    resources:
        mem_mb=4 * 1024,
        java_mem_gb=4 * 1024,
        time_min=get_6h_per_attempt,
        tmpdir=tmp,
    log:
        "logs/gatk/mutect2/{sample}_{version}_{manip}.log",
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