#!/usr/bin/python3.8
# conding: utf-8

"""
Filter a VCF file based on the variant density within a genomic window
"""

from typing import Generator, List, TextIO


def save_lines(lines: List[str],
               density: int,
               out_buffer: TextIO) -> Generator[List[str], None, None]:
    """
    This function edits VCF lines
    """
    for line in lines:
        chomp = line.split("\t")
        if chomp[7] != ".":
            chomp[7] += f";VD={density}"
        else:
            chomp[7] = f"VD={density}"
        out_buffer.write("\t".join(chomp))


# IO paths
vcf_in_path = snakemake.input["vcf"]
vcf_out_path = snakemake.output["vcf"]

# The window size
ws = snakemake.params.get("window_size", 150)
filter_out = snakemake.params.get("filter_vcf", False)

# VCF modifications
update_header = snakemake.params.get("update_header", True)
update_info = snakemake.params.get("update_info", True)
header = f"""##Filtered by VariantDensity.py: ws={ws}
##INFO=<ID=VD,Number=1,Type=Float,Description="Variant density in {ws}bp">\n"""


# Main process
with open(vcf_in_path, "r") as vcfin, open(vcf_out_path, "w") as vcfout:

    # Initialize shared variables
    last_chr = None  # The chromosome that is currently being parsed
    window_init = 0  # The begining of the scanned window
    density = 1  # The number of variant within the scanned window
    lines = []

    for line in vcfin:
        if line.startswith("#"):
            # Then we are reading a header and we shall print it
            if not line.startswith("##") and update_header is True:
                # Then it is the last line of the header
                vcfout.write(header)
            vcfout.write(line)
            continue

        current_chr, pos, *_ = line.split("\t")
        pos = int(pos)

        window_size = pos - window_init <= ws
        same_chr = last_chr == current_chr

        if same_chr and window_size:
            lines.append(line)
            density += 1
            continue

        if (filter_out and density > 1) or not filter_out:
            save_lines(lines, density, vcfout)
        lines = [line]
        last_chr = current_chr
        density = 1
        window_init = pos

    if (filter_out and density > 1) or not filter_out:
        save_lines(lines, density, vcfout)
