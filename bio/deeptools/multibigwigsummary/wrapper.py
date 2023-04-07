__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2023, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

from snakemake.shell import shell
from tempfile import TemporaryDirectory

log = snakemake.log_fmt_shell(stdout=True, stderr=True)
extra = snakemake.params.get("extra", "")


blacklist = snakemake.input.get("blacklist", "")
if blacklist != "":
    blacklist = " --blackListFileName " + blacklist


raw_counts = snakemake.output.get("raw_counts")
if raw_counts != "":
    raw_counts = " --outRawCounts " + raw_counts


subcommand = " bins "
bed = snakemake.input.get("bed", "")
if bed != "":
    bed = " --BED " + bed
    subcommand = " BED-file "


with TemporaryDirectory() as tempdir:
    if "deepBlueURL" in extra:
        extra += "--deepBlueKeepTemp " + tempdir

    shell(
        "multiBigwigSummary {subcommand} "
        "{extra} {blacklist} {raw_counts} "
        "--numberOfProcessors {snakemake.threads} "
        "--bwfiles {snakemake.input.bigwigs} "
        "--outFileName {snakemake.output.npz} "
        "{log} "
    )
