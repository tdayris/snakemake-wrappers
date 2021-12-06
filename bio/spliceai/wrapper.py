# -*- coding: utf-8 -*-

"""Splice AI wrapper"""

from snakemake.shell import shell
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

shell("export LD_LIBRARY_PATH=/usr/local/cuda-11.1/targets/x86_64-linux/lib:${{LD_LIBRARY_PATH}}")
shell("export PATH=/usr/local/cuda/bin:${{PATH}}")
shell("export OMP_NUM_THREADS={snakemake.threads}")

if snakemake.params.piped is True:
    shell(
        "spliceai -R {snakemake.input.fasta} "
        "-A {snakemake.params.annotation} "
        "-O {snakemake.output.vcf} "
        "< {snakemake.input.vcf} {log}"
    )
else:
    shell(
        "spliceai -I {snakemake.input.vcf} "
        "-R {snakemake.input.fasta} "
        "-A {snakemake.params.annotation} "
        "-O {snakemake.output.vcf} {log}"
    )
