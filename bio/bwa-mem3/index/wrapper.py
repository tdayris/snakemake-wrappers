# coding: utf-8

"""Snakemake wrapper for bwa-mem3 index"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2026, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

from os.path import commonprefix, exists
from tempfile import TemporaryDirectory
from snakemake.shell import shell
from shlex import quote
from snakemake_wrapper_utils.snakemake import get_mem

extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)
extra += f" --max-memory '{get_mem(snakemake)}M'"

fasta = str(snakemake.input.fasta)
output_prefix = commonprefix(snakemake.output).rstrip(".")
# Provide output prefix only if needed
if output_prefix != fasta:
    # --meth does not accept -p (outputs <in.fasta>.* and <in.fasta>.meth.*)
    if output_prefix.endswith("meth"):
        output_prefix = output_prefix[:-len(".fasta.meth")]
        extra += " --meth"
        if not exists(f"{output_prefix}.fasta"):
            shell(
                "ln --force --relative --symbolic --verbose "
                "{snakemake.input.fasta:q} {output_prefix}.fasta {log}"
            )
            fasta = f"{output_prefix}.fasta"
    else:
        extra += f" -p {quote(output_prefix)}"
elif output_prefix.endswith("meth"):
    extra += " --meth"

# Ensure compatibility with bwa-mem2 and other tools that require unpacked reference
if any(str(o).endswith("0123") for o in snakemake.output):
    extra += " --emit-unpacked-ref"

with TemporaryDirectory() as tempdir:
    shell(
        "bwa-mem3 index --tmp-dir {tempdir:q} "
        "-t {snakemake.threads} {extra} {fasta:q} {log}"
        " && pwd {log} && tree {log} && ls {log}"
    )
