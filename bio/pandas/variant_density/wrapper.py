#!/usr/bin/python3.8
# conding: utf-8

"""
Filter a VCF file based on the variant density within a genomic window
"""

from typing import Generator, List, TextIO


def save_lines(lines: List[str],
               density: int,
               cluster_id: int,
               window: str,
               out_buffer: TextIO) -> Generator[List[str], None, None]:
    """
    This function edits VCF lines
    """
    for line in lines:
        chomp = line.split("\t")
        if chomp[7] == ".":
            chomp[7] = ""
        else:
            chomp[7] += ";"

        chomp[7] += (
            f"VD={density};"
            f"ClusterID={cluster_id};"
            f"ClusterPosition={window}"
        )
        out_buffer.write("\t".join(chomp))


# IO paths
vcf_in_path = snakemake.input["vcf"]
vcf_out_path = snakemake.output["vcf"]

# The window size
window_size = snakemake.params.get("window_size", 150)
filter_out = snakemake.params.get("filter_vcf", False)

# VCF modifications
update_header = snakemake.params.get("update_header", True)
update_info = snakemake.params.get("update_info", True)
header = f"""##Filtered by VariantDensity.py: ws={window_size}
##INFO=<ID=VD,Number=1,Type=Float,Description="Variant density in {window_size}bp">
##INFO=<ID=ClusterID,Number=1,Type=Float,Description="The name of the variant density cluster">
##INFO=<ID=ClusterPosition,Number=1,Type=String,Description="The cluster position as chr:start-stop">\n"""


# Main process
with open(vcf_in_path, "r") as vcfin, open(vcf_out_path, "w") as vcfout:
    #vcf_iter = iter(vcfin, None)
    window = {
        "chrom": None,
        "start": None,
        "stop": None,
        "density": 0,
        "id": 0
    }
    lines = []

    for line in vcfin:
        if line.startswith("#"):
            # Then we are reading a header and we shall print it
            if not line.startswith("##") and update_header is True:
                # Then it is the last line of the header
                vcfout.write(header)
            vcfout.write(line)
            continue

        chrom, pos, *_ = line.split("\t")
        pos = int(pos)

        if chrom == window["chrom"] and pos - window["start"] <= window_size:
            window["stop"] = pos
            try:
                window['density'] += 1
            except AttributeError:
                window["density"] = 1
            lines.append(line)

        else:
            if (len(lines) != 0) and (len(lines) > 1 or filter_out is False):
                save_lines(
                    lines=lines,
                    density=1 if window["density"] is None else  window["density"],
                    cluster_id=window["id"],
                    window=f"{window['chrom']}:{window['start']}-{window['stop']}",
                    out_buffer=vcfout
                )
            window["chrom"] = chrom
            window["start"] = pos
            window["stop"] = pos
            window["density"] = 1
            window["id"] += 1
            lines = [line]

    # Do not forget the last line!
    if len(lines) > 1 or filter_out is False:
        save_lines(
            lines=lines,
            density=1 if window["density"] is None else  window["density"],
            cluster_id=window["id"],
            window=f"{window['chrom']}:{window['start']}-{window['stop']}",
            out_buffer=vcfout
        )