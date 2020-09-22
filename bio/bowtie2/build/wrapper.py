"""Snakemake wrapper for bowtie2 build"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2020"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"


from snakemake.shell import shell
from os.path import splitext

extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

input = ""
if "fasta" in snakemake.input.keys():
    input = "-f {}".format(snakemake.input["fasta"])
elif "fasta" in snakemake.params.keys():
    input = "-c {}".format(snakemake.params["fasta"])
else:
    raise ValueError(
        "Input sequence could not be found."
    )

prefix = "bwt2_index"
if "prefix" in snakemake.params.keys():
    prefix = snakemake.params["prefix"]


shell(
    " bowtie2-build "
    " {input} "
    " {prefix} "
    " --threads {snakemake.threads} "
    " {extra} "
    " {log} "
)
