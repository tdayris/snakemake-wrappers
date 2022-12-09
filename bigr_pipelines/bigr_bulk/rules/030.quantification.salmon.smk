rule salmon_quantification:
    input:
        **unpack(get_salmon_index),
        **unpack(get_salmon_fastq),
    output:
        quant="030.quantification/salmon_{genome_build}.{genome_release}/{sample}/quant.sf"
    threads: 20
    resources:
        mem_mb=get_10gb_and_10gb_per_attempt,
        time_min=get_45min_and_20min_per_attempt,
        tmpdir="tmp"
    log:
        "logs/030.quantification.salmon/{sample}.{genome_build}.{genome_release}.smk"
    params:
        extra="--numBootstraps 100 --gcBias --seqBias --posBias"
    wrapper:
        "v1.20.0/bio/salmon/quant"