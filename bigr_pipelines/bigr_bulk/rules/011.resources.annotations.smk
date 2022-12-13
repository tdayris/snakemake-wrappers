rule download_GRCh38_108_gtf_annotation:
    output:
        "resources/GRCh38.108.gtf",
    threads: 2
    resources:
        mem_mb=get_512mo_per_attempt,
        time_min=get_2h_per_attempt,
        tmpdir="tmp",
    params:
        species="homo_sapiens",
        datatype="dna",
        build="GRCh38",
        release="108",
    log:
        "logs/011.resources.annotation/GRCh38.108.gtf.log"
    cache: "omit-software"
    wrapper:
        "v1.20.0/bio/reference/ensembl-annotation"


rule download_GRCm39_108_gtf_annotation:
    output:
        "resources/GRCm39.108.gtf",
    threads: 2
    resources:
        mem_mb=get_512mo_per_attempt,
        time_min=get_2h_per_attempt,
        tmpdir="tmp",
    params:
        species="homo_sapiens",
        datatype="dna",
        build="GRCm39",
        release="108",
    log:
        "logs/011.resources.annotation/GRCm39.108.gtf.log"
    cache: "omit-software"
    wrapper:
        "v1.20.0/bio/reference/ensembl-annotation"


rule get_somatic_variants:
    input:
        fai="resources/GRCh38.108.dna.fasta.fai"
    output:
        vcf="resources/GRCh38.108.vcf.gz",
    threads: 1
    resources:
        mem_mb=get_512mo_per_attempt,
        time_min=get_2h_per_attempt,
        tmpdir="tmp",
    params:
        species="homo_sapiens",
        datatype="dna",
        build="GRCh38",
        release="108",
        type="",
    log:
        "logs/011.resources.annotation/GRCh38.108.somatic.log"
    cache: "omit-software"
    wrapper:
        "v1.20.0/bio/reference/ensembl-variation"


rule get_somatic_variants:
    input:
        fai="resources/GRCm39.108.dna.fasta.fai"
    output:
        vcf="resources/GRCm39.108.vcf.gz",
    threads: 1
    resources:
        mem_mb=get_512mo_per_attempt,
        time_min=get_2h_per_attempt,
        tmpdir="tmp",
    params:
        species="homo_sapiens",
        datatype="dna",
        build="GRCm39",
        release="108",
        type="",
    log:
        "logs/011.resources.annotation/GRCm39.108.somatic.log"
    cache: "omit-software"
    wrapper:
        "v1.20.0/bio/reference/ensembl-variation"
