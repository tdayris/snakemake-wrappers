rule bwa_memx_test:
    input:
        "mapped/mem/a.cram",
        "mapped/mem2/a.cram",
        "mapped/meme/a.cram",


rule L2_PARAMETERS_extract:
    input:
        "genome.fasta.suffixarray_uint64_L2_PARAMETERS.gz",
    output:
        temp("genome.fasta.suffixarray_uint64_L2_PARAMETERS"),
    log:
        "log/extract.l2.log",
    conda:
        "gzip.yaml"
    shell:
        "zcat {input} > {output} 2> {log}"
