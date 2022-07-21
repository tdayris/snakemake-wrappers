module snpsift:
    snakefile: str(worflow_source_dir / ".." / ".." / "meta" / "bio" / "snpsift" / "test" / "Snakefile")
    config: snpeff_snpsift_config

use rule * from snpsift