__author__ = "Johannes Köster"
__copyright__ = "Copyright 2016, Johannes Köster"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"


from snakemake.shell import shell
from snakemake_wrapper_utils.samtools import get_samtools_opts


samtools_opts = get_samtools_opts(snakemake)
log = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)

ref = ""
if "fasta" in snakemake.input.keys():
    ref = f"--reference {snakemake.input['fasta']}"

shell(
    "samtools view {snakemake.params.extra} {ref} {samtools_opts} "
    "-o {snakemake.output[0]} {snakemake.input['aln']} "
    "{snakemake.params.position} {log}"
)
