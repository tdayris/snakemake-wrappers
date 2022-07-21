varscan2_somatic_config = {
    "genome": config["ref"]["fasta"],
    "bed": config["ref"]["capture_kit_bed"]
}

module varscan2_somatic_meta:
    snakefile: str(worflow_source_dir / ".." / ".." / "meta" / "bio" / "varscan2_somatic" / "test" / "Snakefile")
    config: varscan2_somatic_config

use rule * from varscan2_somatic_meta