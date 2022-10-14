__author__ = "Jan Forster"
__copyright__ = "Copyright 2020, Jan Forster"
__email__ = "j.forster@dkfz.de"
__license__ = "MIT"


import tempfile
from pathlib import Path
from snakemake.shell import shell
from snakemake_wrapper_utils.bcftools import get_bcftools_opts


bcftools_opts = get_bcftools_opts(snakemake, parse_ref=False, parse_memory=False)
extra = snakemake.params.get("extra", "")
view_extra = snakemake.params.get("view_extra", "")
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

if snakemake.threads != 2:
    raise ValueError("BCFTools reheader wrapper requires exactly 2 threads")

## Extract arguments
header = snakemake.input.get("header", "")
if header:
    header = f"-h {header}"

samples = snakemake.input.get("samples", "")
if samples:
    samples = f"-s {samples}"

fai = snakemake.input.get("fai", None)
if fai is not None:
    fai_cmd = f"--fai {fai}"
else:
    fai_cmd = ""


with tempfile.TemporaryDirectory() as tmpdir:
    tmp_prefix = Path(tmpdir) / "bcftools_reheader."


shell(
    "bcftools reheader "
    "{extra} "
    "{header_cmd} "
    "{samples_cmd} "
    "{fai_cmd} "
    "{snakemake.input.vcf} "
    "| bcftools view "
    "{view_extra} "
    " {bcftools_opts}"
    "> {snakemake.output} "
    "{log}"
)
