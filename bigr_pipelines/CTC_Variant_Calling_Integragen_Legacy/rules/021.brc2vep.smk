"""
Rscript provided
"""


rule brc2vep_ctc:
    input:
        readcount="bam_readcount/{sample}_{version}_{manip}_{nb}.tsv",
    output:
        tumor_dp20=temp("bcr2vep/dp/{sample}_{version}_{manip}_{nb}.tsv"),
        filtered=temp("bcr2vep/filtered/{sample}_{version}_{manip}_{nb}.tsv"),
    threads: 1
    resources:
        mem_mb=get_10gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir=tmp,
    log:
        "logs/bcr2vep/{sample}_{version}_{manip}_{nb}.log",
    params:
        dp=config["params"].get("min_dp", 20),
        alt_ad=config["params"].get("min_alt_ad", 10),
        alt_vaf=config["params"].get("min_alt_vaf", 0.03),
    conda:
        str(workflow_source_dir / "envs" / "r.yaml")
    script:
        str(workflow_source_dir / "scripts" / "bcr2vep.R")



rule brc2vep_wbc:
    input:
        readcount="bam_readcount/{sample}_{version}_{manip}.tsv",
    output:
        tumor_dp20=temp("bcr2vep/dp/{sample}_{version}_{manip}.tsv"),
        filtered=temp("bcr2vep/filtered/{sample}_{version}_{manip}.tsv"),
    threads: 1
    resources:
        mem_mb=get_10gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir=tmp,
    log:
        "logs/bcr2vep/{sample}_{version}_{manip}.log",
    params:
        dp=config["params"].get("min_dp", 20),
        alt_ad=config["params"].get("min_alt_ad", 10),
        alt_vaf=config["params"].get("min_alt_vaf", 0.03),
    conda:
        str(workflow_source_dir / "envs" / "r.yaml")
    script:
        str(workflow_source_dir / "scripts" / "bcr2vep.R")




rule brc2vep_baseline:
    input:
        readcount="bam_readcount/{sample}.baseline.tsv",
    output:
        tumor_dp20=temp("bcr2vep/dp/{sample}.baseline.tsv"),
        filtered=temp("bcr2vep/filtered/{sample}.baseline.tsv"),
    threads: 1
    resources:
        mem_mb=get_10gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir=tmp,
    log:
        "logs/bcr2vep/{sample}.baseline.log",
    params:
        dp=config["params"].get("min_dp", 20),
        alt_ad=config["params"].get("min_alt_ad", 10),
        alt_vaf=config["params"].get("min_alt_vaf", 0.03),
    conda:
        str(workflow_source_dir / "envs" / "r.yaml")
    script:
        str(workflow_source_dir / "scripts" / "bcr2vep.R")
