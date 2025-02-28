
"""
No command line provided.
"""

rule zip_variants_hc_baseline:
    input:
        "gatk/select_variants/baseline/{sample}.vcf"
    output:
        protected("data_output/HC_Baseline/{sample}.vcf.gz")
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


rule zip_variants_hc_wbc:
    input:
        "gatk/select_variants/wbc/{sample}_{version}_{manip}.vcf"
    output:
        protected("data_output/HC_WBC/{sample}_{version}_{manip}.vcf.gz")
    threads: 2
    resources:
        mem_mb=get_4gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir="tmp"
    log:
        "logs/bcftools/view/{sample}_{version}_{manip}.wbc.log"
    params:
        extra=""
    wrapper:
        "bio/bcftools/view"



rule zip_variants_hc_ctc:
    input:
        "gatk/select_variants/ctc/{sample}_{version}_{manip}_{nb}.vcf"
    output:
        protected("data_output/HC_CTC/{sample}_{version}_{manip}_{nb}.vcf.gz")
    threads: 2
    resources:
        mem_mb=get_4gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir="tmp"
    log:
        "logs/bcftools/view/{sample}_{version}_{manip}_{nb}.ctc.log"
    params:
        extra=""
    wrapper:
        "bio/bcftools/view"


rule tabix_variants_hc_baseline:
    input:
        "data_output/HC_Baseline/{sample}.vcf.gz"
    output:
        protected("data_output/HC_Baseline/{sample}.vcf.gz.tbi")
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


rule tabix_variants_hc_wbc:
    input:
        "data_output/HC_WBC/{sample}_{version}_{manip}.vcf.gz"
    output:
        protected("data_output/HC_WBC/{sample}_{version}_{manip}.vcf.gz.tbi")
    threads: 1
    resources:
        mem_mb=get_4gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir="tmp"
    log:
        "logs/tabix/{sample}_{version}_{manip}.wbc.log"
    params:
        "-p vcf"
    wrapper:
        "bio/tabix/index"


rule tabix_variants_hc_ctc:
    input:
        "data_output/HC_CTC/{sample}_{version}_{manip}_{nb}.vcf.gz"
    output:
        protected("data_output/HC_CTC/{sample}_{version}_{manip}_{nb}.vcf.gz.tbi")
    threads: 1
    resources:
        mem_mb=get_4gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir="tmp"
    log:
        "logs/tabix/{sample}_{version}_{manip}_{nb}.ctc.log"
    params:
        "-p vcf"
    wrapper:
        "bio/tabix/index"