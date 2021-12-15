
rule splice_ai:
    input:
        vcf = "snpsift/clinvar/{sample}.vcf",
        fasta = "/mnt/beegfs/database/bioinfo/Index_DB/Fasta/Ensembl/GRCh38.99/GRCh38.99.homo_sapiens.dna.fasta"
    output:
        vcf = temp("splice_ai/{sample}.vcf")
    message:
        "Adding Splice Variant annotation to {wildcards.sample}"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 2048 * 5,
        time_min=lambda wildcards, attempt: attempt * 60 * 1.5,
        tmpdir="tmp",
        gres="gpu:t4:1"
    params:
        annotation = "grch38",
        piped = True
    log:
        "logs/splice_ai/{sample}.log"
    wrapper:
        "bio/spliceai"
