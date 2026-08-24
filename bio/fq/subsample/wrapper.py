# coding: utf-8

"""snakemake wrapper for fq lint"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2026, Thibault Dayris"
__license__ = "MIT"

from snakemake.shell import shell
from shlex import quote

log = snakemake.log_fmt_shell(stdout=True, stderr=True)
extra = snakemake.params.get("extra", "")

fq_nb = len(snakemake.input)
if not (1 <= fq_nb <= 2):
    raise ValueError(f"This wrapper expected 1 or 2 fastq files, got {fq_nb}.")

out_nb = len(snakemake.output)
if fq_nb != out_nb:
    raise ValueError(
        "This wrapper expected as many input/output files. "
        f"Got {fq_nb} input, and {out_nb} output."
    )

extra += f" --r1-dst {quote(snakemake.output[0])}"
if out_nb == 2:
    extra += f" --r2-dst {quote(snakemake.output[1])}"

shell("fq subsample {extra} {snakemake.input} {log}")

