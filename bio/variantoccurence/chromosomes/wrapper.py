"""Snakemake wrapper to extract variants and count occurences"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2021, Dayris Thibault"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"


from snakemake.shell import shell
import os.path as op

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

if snakemake.threads != 7:
    raise ValueError("7 threads are required for this wrapper.")

if all(i.endswith("vcf") for i in snakemake.input["calls"]):
    reader = "cat"
elif all(i.endswith("vcf.gz") for i in snakemake.input["calls"]):
    reader = "zcat"
else:
    raise ValueError(
        "All VCF should be bgzipped, or all VCF should be raw text"
    )

chr = snakemake.wildcards["chr"]

shell(
    """(echo -e "Occurence\tChromosome\tPosition\tID\tRef\tAlt"; """
    """for CALL in {snakemake.input.calls} ; do """
    """{reader} ${{CALL}} | cut -f -5 | grep -P "^{chr}" | """
    """ sort | uniq ; done | sort | uniq -c ) """
    """> {snakemake.output.txt} {log}"""
)
