"""
No command line provided
"""


rule brc2vep:
    input:
        readcount="bam_readcount/{sample}.{status}.tsv",
    output:
        dp20=temp("bcr2vep/dp/{sample}.csv"),
        filtered=temp("bcr2vep/filtered/{sample}.csv"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 10,
        time_min=lambda wildcards, attempt: attempt * 45,
        tmpdir="tmp",
    log:
        "logs/bcr2vep/{sample}.log",
    params:
        dp=config["params"].get("min_dp", 20),
        alt_ad=config["params"].get("min_alt_ad", 10),
        alt_vaf=config["params"].get("min_alt_vaf", 0.03),
    conda:
        "envs/r.yaml"
    script:
        "scripts/bcr2vep.R"
