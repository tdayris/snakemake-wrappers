#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Format Revel CSV to VCF"""

import datetime
import logging

from typing import Optional

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

def build_info(**kwargs) -> str:
    return ";".join([f"REVEL_{k}={v}" for k, v in kwargs.items()])

def vcf_entry(
        chrom: str,
        pos: int,
        vatid: Optional[str] = ".",
        ref: Optional[str] = ".",
        alt: Optional[str] = ".",
        qual: Optional[str] = ".",
        filter: Optional[str] = ".",
        info: Optional[str] = ".",
        fmt: Optional[str] = None,
        *samples: Optional[str] = None
    ) -> str:
    """Create a VCF entry line"""
    if fmt is None:
        return "\t".join([chrom, pos, varid, ref, alt, qual, filter, info]) + "\n"
    return "\t".join([chrom, pos, varid, ref, alt, qual, filter, info, fmt, *samples]) + "\n"


version = 1.0
name = "Revel_to_vcf"
url = "github.com/tdayris/snakemake-wrappers/tree/Unofficial/bio/BiGR/revel/tovcf/wrapper.py"
headers = [
    f"""##BiGRCommandLine=<ID={name},CommandLine="{url}",Version={version},Date={datetime.date.today()}>\n""",
    f"""##INFO=<ID=REVEL_chr,Number=1,Type=String,Description="Revel chromosome ID">\n""",
    f"""##INFO=<ID=REVEL_hg19_pos,Number=1,Type=Integer,Description="Revel HG19 Position">\n""",
    f"""##INFO=<ID=REVEL_grch38_pos,Number=1,Type=Integer,Description="Revel HG38 Position">\n""",
    f"""##INFO=<ID=REVEL_ref,Number=1,Type=String,Description="Revel reference base">\n""",
    f"""##INFO=<ID=REVEL_alt,Number=1,Type=String,Description="Revel alternative site">\n""",
    f"""##INFO=<ID=REVEL_aaref,Number=1,Type=String,Description="Revel reference aminoacid">\n""",
    f"""##INFO=<ID=REVEL_aaalt,Number=1,Type=String,Description="Revel alternative aminoacid">\n""",
    f"""##INFO=<ID=REVEL_REVEL,Number=1,Type=Float,Description="Revel score">\n""",
    f"""##INFO=<ID=REVEL_Ensembl_transcriptid,Number=1,Type=String,Description="Revel ensembl transcript ID">\n"""
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n
]


if str(snakemake.output["vcf"]).endswith("vcf"):
    output_file = snakemake.output["vcf"]
else:
    output_file = snakemake.output["vcf"][:-3]

with (open(snakemake.input["db"], "r") as revel,
      open(output_file, "w") as vcf):
    vcf.write("".join(headers))
    fields = next(revel)[:-1].split(",")
    for line in revel:
        chomp=line[:-1].split(",")
        entry = vcf_entry(
            chrom=chomp[0], pos=chomp[2],
            ref=chomp[3], alt=chomp[4],
            info=build_info(**dict(zip(fields, chomp))
        )
        vcf.write(entry)

if str(snakemake.output["vcf"]).endswith(".gz"):
    shell("pbgzip -c {uncompressed_vcf} > {snakemake.output['vcf']} 2> {log}")
    shell("tabix -p vcf {snakemake.output['vcf']} >> {log} 2>&1")
    shell("rm --verbose uncompressed_vcf >> {log} 2>&1")
