import pandas

design = pandas.read_table(
    config["design"], sep="\t", header=0,
).set_index(
    keys="Sample_id", drop=False, verify_integrity=True
)

if not "Protocol" in snakemake.columns.tolist():
    logging.warning(
        "Protocol was not set in design, "
        "all samples are treated as rna-seq bulk samples"
    )


if not "Capturekit_bed" in snakemake.columns.tolist():
    logging.warning(
        "Capturekit_bed was not set in design, "
        "all samples are treated as belonging to "
        "the same (default) capturekit bed."
    )


if not "Organism" in snakemake.columns.tolist():
    logging.warning(
        "Organism was not set in design, "
        "all samples are treated as belonging to "
        "the same organism."
    )