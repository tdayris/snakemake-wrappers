from snakemake.utils import validate

configfile: str(workflow_source_dir / "config" / "config.yaml")
validate(config, str(workflow_source_dir / "schema" / "config.schema.yaml"))