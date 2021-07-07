#!/usr/bin/env python

"""
Snakemake wrapper for cnv_facets
"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2020, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"


from snakemake.shell import shell
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

extra = params.get("extra", "")
if "bed" in snakemake.input:
    extra += " --targets {input.bed}"

prefix = snakemake.output["vcf"][:-len(".vcf.gz")]

if ("--no-cov-plot" in extra) and ("coverage" in snakemake.output.keys()):
    raise ValueError(
        "Coverage plot was deactivated in extra parameters, but "
        "expected as result. Remove the extra parameter or the "
        "coverage file in output."
    )

if all(i in snakemake.input.keys() for i in ["tumor_bam", "normal_bam", "vcf"]):
    shell(
        "cnv_facets.R "
        "--snp-tumour {input.tumor_bam} "
        "--snp-normal {input.normal_bam} "
        "--snp-vcf {input.vcf} "
        "--out {prefix} "
        "--snp-nprocs {threads} "
        "{extra} "
        "{log}"
    )
elif "pileup" in snakemake.input.keys():
    shell(
        "cnv_facets.R "
        "--pileup {input.pileup} "
        "--out {output.prefix} "
        "{extra} "
        "{log}"
    )
else:
    raise KeyError(
        "Expecting either: a tumor and normal bam with reference VCF, "
        "or a pileup file."
    )
