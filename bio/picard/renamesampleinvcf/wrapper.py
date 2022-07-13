"""Snakemake wrapper for picard RenameSampleInVcf."""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2020, Dayris Thibault"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"


from snakemake.shell import shell
from snakemake_wrapper_utils.java import get_java_opts

log = snakemake.log_fmt_shell(stdout=False, stderr=True)
extra = snakemake.params.get("extra", "")
java_opts = get_java_opts(snakemake)

if "index" in snakemake.output.keys():
    extra += " --CREATE_INDEX true"

if "fasta" in snakemake.input.keys():
    extra += "--REFERENCE_SEQUENCE {}".format(snakemake.input["fasta"])


shell(
    "picard"
    " RenameSampleInVcf"
    " {java_opts}"
    " {extra}"
    " INPUT={snakemake.input.vcf}"
    " OUTPUT={snakemake.output.vcf}"
    " NEW_SAMPLE_NAME={snakemake.params.new_sample_name}"
    " {log}"
)
