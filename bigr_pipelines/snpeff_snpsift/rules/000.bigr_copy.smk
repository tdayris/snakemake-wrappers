# Copy files present in the design file
# takes iRODS input as well as absolute paths
rule bigr_copy:
    output:
        "data_input/calls/{sample}.vcf.gz"
    message:
        "Getting {wildcards.sample} VCF file"
    threads: 1
    resources:
      mem_mb=lambda wildcards, attempt: min(attempt * 1024, 2048),
      time_min=lambda wildcards, attempt: attempt * 45,
    params:
        input=lambda wildcards, output: vcf_links[output[0]]
    log:
        "logs/bigr_copy/{sample}.log"
    wrapper:
        "bio/BiGR/copy"