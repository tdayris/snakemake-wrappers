rule splice_ai:
    input:
        vcf="bcftools/filter/{sample}.preannot.vcf",
        fasta=config["reference"]["fasta"],
        fasta_index=config["reference"]["fasta_index"],
    output:
        vcf=temp("splice_ai/annot/{sample}.vcf.gz"),
    threads: config.get("max_threads", 20)
    resources:
        mem_mb=get_8gb_per_attempt,
        time_min=get_6h_per_attempt,
        tmpdir=tmp,
    retries: 1
    params:
        annotation=config["reference"].get("ncbi_build", "grch38").lower(),
        piped=True,
    log:
        "logs/splice_ai/{sample}.log",
    wrapper:
        str(wrapper_prefix / "bio" / "spliceai")
