rule salmon_build_gentrome:
    input:
        **unpack(get_gentrome),
    output:
        gentrome=temp("resources/gentrome.fasta"),
        decoys=temp("resources/decoys.txt"),
    threads: 1
    resources:
        mem_mb=get_512mo_per_attempt,
        time_min=get_10min_per_go,
        tmpdir="tmp",
    log:
        "logs/012.resources/salmon_build_gentrome.smk"
    wrapper:
        "v1.20.0/bio/salmon/decoys"


rule salmon_index_decoy_aware_transcriptome:
    input:
        sequences="resources/{genome_build}.{genome_release}.gentrome.fasta",
        decoy="resources/{genome_build}.{genome_release}.decoys.txt",
    output:
        multiext(
            "resources/salmon/{genome_build}.{genome_release}.index/",
            "complete_ref_lens.bin",
            "ctable.bin",
            "ctg_offsets.bin",
            "duplicate_clusters.tsv",
            "info.json",
            "mphf.bin",
            "pos.bin",
            "pre_indexing.log",
            "rank.bin",
            "refAccumLengths.bin",
            "ref_indexing.log",
            "reflengths.bin",
            "refseq.bin",
            "seq.bin",
            "versionInfo.json",
        ),
    threads: 20
    resources:
        mem_mb=get_10gb_and_10gb_per_attempt,
        time_min=get_1h_per_attempt,
        tmpdir="tmp",
    params:
        extra="",
    log:
        "logs/012.resources/salmon_index_decoy_aware_transcriptome/{genome_build}.{genome_release}.log",
    cache: True
    wrapper:
        "v1.20.0/bio/salmon/index"
