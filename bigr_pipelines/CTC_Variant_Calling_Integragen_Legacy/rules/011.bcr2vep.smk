"""
No command line provided
"""


rule brc2vep:
    input:
        readcount="bam_readcount/{sample}.ctc.tsv",
    output:
        tumor_dp20=temp("bcr2vep/dp/{sample}.ctc.tsv"),
        filtered=temp("bcr2vep/filtered/{sample}.ctc.tsv"),
    threads: 1
    resources:
        mem_mb=get_10gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir="tmp",
    log:
        "logs/bcr2vep/{sample}.log",
    params:
        dp=config["params"].get("min_dp", 20),
        alt_ad=config["params"].get("min_alt_ad", 10),
        alt_vaf=config["params"].get("min_alt_vaf", 0.03),
    conda:
        str(workflow_source_dir / "envs" / "r.yaml")
    script:
        "scripts/bcr2vep.R"
