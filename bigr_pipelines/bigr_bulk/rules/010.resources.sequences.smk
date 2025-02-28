rule download_GRCh38_108_dna_sequence:
    output:
        "resources/GRCh38.108.dna.fasta",
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
        "logs/010.resources.sequences/GRCh38.108.dna.log"
    cache: "omit-software"
    wrapper:
        "v1.20.0/bio/reference/ensembl-sequence"


rule download_GRCm39_108_dna_sequence:
    output:
        "resources/GRCm39.108.dna.fasta",
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
        "logs/010.resources.sequences/GRCm39.108.dna.log"
    cache: "omit-software"
    wrapper:
        "v1.20.0/bio/reference/ensembl-sequence"


rule download_GRCh38_108_cdna_sequence:
    output:
        "resources/GRCh38.108.cdna.fasta",
    threads: 2
    resources:
        mem_mb=get_512mo_per_attempt,
        time_min=get_2h_per_attempt,
        tmpdir="tmp",
    params:
        species="homo_sapiens",
        datatype="cdna",
        build="GRCh38",
        release="108",
    log:
        "logs/010.resources.sequences/GRCh38.108.cdna.log"
    cache: "omit-software"
    wrapper:
        "v1.20.0/bio/reference/ensembl-sequence"


rule download_GRCm39_108_cdna_sequence:
    output:
        "resources/GRCm39.108.cdna.fasta",
    threads: 2
    resources:
        mem_mb=get_512mo_per_attempt,
        time_min=get_2h_per_attempt,
        tmpdir="tmp",
    params:
        species="homo_sapiens",
        datatype="cdna",
        build="GRCm39",
        release="108",
    log:
        "logs/010.resources.sequences/GRCm39.108.cdna.log"
    cache: "omit-software"
    wrapper:
        "v1.20.0/bio/reference/ensembl-sequence"


rule samtools_faidx_reference:
    input:
        ancient("resources/{genome_build}.{genome_release}.{seqtype}.fasta"),
    output:
        "resources/{genome_build}.{genome_release}.{seqtype}.fasta.fai",
    threads: 1
    resources:
        mem_mb=get_input_size,
        time_min=get_10min_per_go,
        tmpdir="tmp",
    params:
        extra="",
    log:
        "logs/010.resources.sequences/samtools_faidx/{genome_build}.{genome_release}.{seqtype}.log"
    cache: True
    wrapper:
        "v1.20.0/bio/samtools/faidx"
    

rule picard_dict_reference:
    input:
        ancient("resources/{genome_build}.{genome_release}.{seqtype}.fasta"),
    output:
        "resources/{genome_build}.{genome_release}.{seqtype}.dict",
    threads: 1
    resources:
        mem_mb=get_input_size_plus_3go,
        time_min=get_10min_per_go,
        tmpdir="tmp",
    params:
        extra="",
    log:
        "logs/010.resources.sequences/picard_dict_reference/{genome_build}.{genome_release}.{seqtype}.log"
    cache: True
    wrapper:
        "v1.20.0/bio/picard/createsequencedictionary"