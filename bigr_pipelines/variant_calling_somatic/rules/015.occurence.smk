rule variant_occurence_annotate:
    input:
        calls=["snpsift/dbvar/{sample}.vcf"],
        occurence="bigr/occurences/all_chroms.txt",
    output:
        calls=[temp("bigr/occurence_annotated/{sample}.vcf")],
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
        time_min=lambda wildcards, attempt: attempt * 15,
        tmpdir="tmp",
    log:
        "logs/variant_occurence/uncompress/{sample}.log",
    wrapper:
        "bio/variantoccurence/annotate"


rule concatenate_per_chr_information:
    input:
        expand("bigr/occurence/{chr}.txt", chr=config["params"]["chr"]),
    output:
        temp("bigr/occurences/all_chroms.txt"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
        time_min=lambda wildcards, attempt: attempt * 15,
        tmpdir="tmp",
    log:
        "logs/variant_occurence/all.log",
    shell:
        "for i in {input}; do sed '1d' ${{i}}; done > {output} 2> {log}"


rule variant_occurence_per_chr:
    input:
        calls=expand("snpsift/vartype/{sample}.vcf", sample=samples_list),
    output:
        txt=temp("bigr/occurence/{chr}.txt"),
    threads: 7
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
        time_min=lambda wildcards, attempt: attempt * 45,
        tmpdir="tmp",
    log:
        "logs/variant_occurence/{chr}.log",
    wrapper:
        "bio/variantoccurence/chromosomes"
