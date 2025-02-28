rule subset_gene_counts:
    input:
        table="salmon/TPM.genes.tsv"
    output:
        table="immunedeconv/TPM.tsv"
    message:
        "Formatting counts for ImmuneDeconv"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 512,
        time_min=lambda wildcards, attempt: attempt * 15,
        tmpdir="tmp"
    log:
        "logs/immunedeconv/filter_gene_counts.log"
    params:
        drop_column = ["target_id"],
        set_index = "Hugo_ID",
        drop_duplicated_lines = True,
        keep_index = True,
        drop_duplicated_index = True,
        override_previous_index = True,
        filters = [
            ["Hugo_ID", "!=", "Unknown"]
        ]
    wrapper:
        "bio/pandas/filter_table"