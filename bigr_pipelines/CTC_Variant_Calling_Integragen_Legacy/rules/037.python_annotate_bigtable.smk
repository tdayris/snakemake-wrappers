"""
No command line provided
"""

rule bigtable_annotated:
    input:
        bigtable="bigtable/bigtable.uniq.csv",
        label="Labels.csv"
    output:
        bigtable="data_output/BigTable.tsv",
        fulltable="upload/BigTable.tsv",
        excel="upload/BigTable.xlsx"
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_15min_per_attempt,
        tmpdir=tmp,
    log:
        "logs/bigtable/annot.log",
    params:
        "",
    conda:
        str(workflow_source_dir / "envs" / "bash.yaml")
    script:
        str(workflow_source_dir / "scripts" / "annotate_bigtable.py")