__author__ = "Johannes Köster"
__copyright__ = "Copyright 2016, Johannes Köster"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"


from snakemake.shell import shell
from snakemake_wrapper_utils.samtools import get_samtools_opts


samtools_opts = get_samtools_opts(snakemake)
log = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)
aln = snakemake.input.get("aln", None)
if aln is None:
    aln = snakemake.input[0]

ref = ""
if "fasta" in snakemake.input.keys():
    ref = f"--reference {snakemake.input['fasta']}"

if "fasta_idx" in snakemake.input.keys():
    ref += f" --reference {snakemake.input['fasta_idx']}"

bed = ""
if "bed" in snakemake.input.keys():
    bed = f"-L {snakemake.input['bed']}"
    
position = snakemake.params.get("position", "")

shell(
    "samtools view {snakemake.params.extra} {ref} {bed} {samtools_opts} "
    "-o {snakemake.output[0]} {aln} "
    "{position} {log}"
)
