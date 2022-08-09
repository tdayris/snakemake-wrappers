__author__ = "Johannes Köster"
__copyright__ = "Copyright 2016, Johannes Köster"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"


import os

from snakemake.shell import shell
from snakemake_wrapper_utils.samtools import get_samtools_opts


samtools_opts = get_samtools_opts(snakemake)
log = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)
aln = snakemake.input.get("aln", None)
if aln is None:
    aln = snakemake.input[0]

ref = snakemake.input.get('fasta')
if ref:
    ref_dir = os.path.dirname(ref)
    ref = f"-T {ref} "

    if "fasta_idx" in snakemake.input.keys():
        ref += f" -t {snakemake.input['fasta_idx']} "

    cache = snakemake.input.get("cache", "")
    if cache:
       os.environ["REF_PATH"] = f"{cache}/%2s/%2s/%s:{ref_dir}:http://www.ebi.ac.uk/ena/cram/md5/%s"
       os.environ["REF_CACHE"] = f"{cache}/%2s/%2s/%s"

bed = ""
if "bed" in snakemake.input.keys():
    bed = f"-L {snakemake.input['bed']} "
    
position = snakemake.params.get("position", "")

shell(
    "samtools view {snakemake.params.extra} {ref} {bed} {samtools_opts} "
    "-o {snakemake.output[0]} {aln} "
    "{position} {log}"
)
