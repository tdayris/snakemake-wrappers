#!/usr/bin/python3.9
# -*- coding: utf-8 -*-

"""
Annotate VCF with fields dedicated to GLeaves
"""

import datetime
import gzip
import logging

from typing import Any, Dict, Optional

logging.basicConfig(
    #filename=snakemake.log[0],
    #filemode="w",
    level=logging.DEBUG
)


def open_function(file: str):
    """Return the correct opening function"""
    if file.endswith(".gz"):
        return gzip.open(file, "rb")
    return open(file, "r")



def get_headers(description: Dict[str, Any]) -> str:
    """
    From a list of column name, and an optional list of description,
    build VCF headers.
    """
    version = 1.0
    name = "gleaves_compatibility"
    url = f"github.com/tdayris/snakemake-wrappers/tree/Unofficial/bio/BiGR/{name}/wrapper.py"

    headers = [
        f"""##BiGRCommandLine=<ID<{name},CommandLine="{url}",Version={version},Date={datetime.date.today()}>\n"""
    ]

    for key, value in description.items():
        headers.append(
            f"""##INFO=<ID={key},Number={value["nb"]},Type={value["type"]},Description="{value["desc"]}">\n"""
        )
    
    return "".join(headers)


headers_description = {
    "VAR": {"nb": "1", "type": "String", "desc": "Genotype Quality (LSC: Likely subclonal (VAF [0-10%]), PSC: Probably subclonal (VAF [10%-30%]; PHE: Probably Heterozygote (VAF [30%-40%] ou [60%-70%]); LHE : Likely Heterozygote (VAF [40%-60%]), PHO: Probably Homozygote (VAF [70%-80%]), LHO: Likely Homozygote (VAF > 80%))"},
    "TRC": {"nb": "1", "type": "Integer", "desc": "Total Read Count"},
    "RRC": {"nb": "1", "type": "Integer", "desc": "Reference allele depth"},
    "ARC": {"nb": "1", "type": "Integer", "desc": "Variant allele depth"},
    "BRC": {"nb": "1", "type": "Integer", "desc": "Number of reads neither on reference, nor on alternative allele"},
    "ARR": {"nb": "1", "type": "Float", "desc": "Allele frquency"},
    "ARCp": {"nb": "1", "type": "Integer", "desc": "ALT supporting depth on positive strand"},
    "ARCm": {"nb": "1", "type": "Integer", "desc": "ALT supporting depth on negative strand"},
    "RRCp": {"nb": "1", "type": "Integer", "desc": "REF supporting depth on positive strand"},
    "RRCm": {"nb": "1", "type": "Integer", "desc": "REF supporting depth on negative strand"},
    "BRR": {"nb": "1", "type": "Float", "desc": "Background noise ratio"},
    "BRE": {"nb": "1", "type": "Float", "desc": "Background noise enrichment"},
    "SBP": {"nb": "1", "type": "Float", "desc": "Background noise estimated using Fisher's exact test"},
    "SBM": {"nb": "1", "type": "Float", "desc": "Background noise estimated by the symmetric odds ratio test"},
    "BKG": {"nb": "1", "type": "String", "desc": "Variant quality (LCL : Likely Clean (BRE [0-20%]); PCL : Probably Clean (BRE [20%-30%]); PNO : Probably Noisy (BRE [30%-50%]); LNO : Likely Noisy (BRE > 80%)"},
    "OLD_MULTIALLELIC": {"nb": "1", "type": "String", "desc": "Original chr:pos:ref:alt encoding"},
    "OLD_VARIANT": {"nb": "1", "type": "String", "desc": "Original chr:pos:ref:alt encoding"}
}


def parse_info(chrom: str, pos: int, ref: str, alt: str, info: str) -> Dict[str, Any]:
    """
    Parse info fields and gather annotations
    """
    new_annotations = {
        "ARCp": None, "ARCm": None, "RRCp": None, "RRCm": None, "VAR": None,
        "SBM": None, "SBP": None, 
        "TRC": None, "BRC": None, "BRR": None, "BRE": None, "BKG": None,
        "OLD_MULTIALLELIC": None,
        "OLD_VARIANT": None
    }
    fields = info.split(";")
    for pos, field in enumerate(fields):
        if field.upper().startswith("AS_SB_TABLE"):
            new_annotations.update(**get_sb_table(field))
        elif field.upper().startswith("SOR"):
            new_annotations.update(**get_SOR(field))
        elif field.upper().startswith("FS"):
            new_annotations.update(**get_fisher_test(field))
        elif field.upper().startswith("AF"):
            new_annotations.update(**get_var(field))
        elif field.upper().startswith("DP"):
            new_annotations.update(**get_trc(field))

    try:
        new_annotations["RRC"] = new_annotations["RRCp"] + new_annotations["RRCm"]
    except TypeError:
        pass


    try:
        new_annotations["ARC"] = new_annotations["ARCp"] + new_annotations["ARCm"]
    except TypeError:
        pass

    try:
        new_annotations["BRC"] = new_annotations["TRC"] - (new_annotations["RRC"] + new_annotations["ARC"])
    except TypeError:
        pass

    try:
        new_annotations["BRR"] = new_annotations["BRC"] / new_annotations["TRC"]
    except TypeError:
        pass

    try:
        new_annotations["BRE"] = 1 - new_annotations["BRR"]
    except TypeError:
        pass


    new_annotations["BKG"] = get_bkg(new_annotations["BRE"])
    new_annotations["OLD_MULTIALLELIC"] = f"{chrom}:{pos}:{ref}:{alt}"
    new_annotations["OLD_VARIANT"] = f"{chrom}:{pos}:{ref}:{alt}"
    return new_annotations
        


def get_sb_table(sb_table: str) -> Dict[str, int]:
    """Parse INFO field and return strand bias table"""
    sb_table = sb_table.split("=")[-1]
    alt, ref = sb_table.split("|")
    ARCp, ARCm = map(int, alt.split(","))
    RRCp, RRCm = map(int, ref.split(","))
    return {"ARCp": ARCp, "ARCm": ARCm, "RRCp": RRCp, "RRCm": RRCm}


def get_SOR(sor_field: str) -> Dict[str, float]:
    """
    Parse INFO field and return Strand Odds Ratio
    """
    return {"SBM": sor_field.split("=")[-1]}


def get_fisher_test(fisher_field: str) -> Dict[str, float]:
    """
    Parse INFO field and return Fisher Exact test
    """
    return {"SBP": fisher_field.split("=")[-1]}


def get_var(vaf_field: str) -> Dict[str, str]:
    """
    Parse VAF field to return genotype quality
    """
    vaf = float(vaf_field.split("=")[-1])
    if 0 <= vaf <= 0.1:
        return {"VAR": "LSC"}
    if 0.1 < vaf < 0.3:
        return {"VAR": "LSC"}
    if (0.3 <= vaf <= 0.4) or (0.6 <= vaf <= 0.7):
        return {"VAR": "PHE"}
    if 0.4 < vaf < 0.6:
        return {"VAR": "LHE"}
    if 0.7 < vaf <= 0.8:
        return {"VAR": "PHO"}
    return {"VAR": "LHO"}


def get_trc(depth: str) -> int:
    """Parse INFO DP and return TRC"""
    return {"TRC": depth.split("=")[-1]}


def get_bkg(bre: Optional[float]) -> Optional[str]:
    """Render BRE human-readably"""
    if bre is None:
        return None

    if 0 <= bre <= 0.2:
        return "LCL"
    if 0.2 < bre <= 0.3:
        return "PCL"
    if 0.3 < bre <= 0.5:
        return "PNO"
    if 0.8 >= bre:
        return "LNO"

    return None


def annotate(line: str) -> str:
    line = line[:-1].split("\t")
    new_annotations = parse_info(
        chrom=line[0],
        pos=line[1],
        ref=line[3],
        alt=line[4],
        info=line[7]
    )
    new_annotations = [
        "{key}={val}" for key, val in new_annotations.items()
        if val is not None
    ]
    new_annotations = ";".join(new_annotations)
    if line[7] == ".":
        if new_annotations != "":
            line[7] = new_annotations
    elif new_annotations != "":
        line[7] += ";"
        line[7] += new_annotations

    if line[6] == ".":
        line[6] = "PASS"
    
    return "\t".join(line) + "\n"

    

# Annotating input VCF
logging.debug("Opening VCF")
if str(snakemake.output["vcf"]).endswith("vcf.gz"):
    out_vcf = snakemake.output["vcf"][:-3]
else:
    out_vcf = snakemake.output["vcf"]



with (open_function(snakemake.input["vcf"]) as in_vcf,
      open(out_vcf, 'w', encoding="utf-8") as out_vcf):
    for line in in_vcf:
        if isinstance(line, bytes):
            line = line.decode("utf-8")

        if line.startswith("#"):
            pass
        else:
            line = annotate(line)

        out_vcf.write(line)


if str(snakemake.output["vcf"]).endswith("vcf.gz"):
    logging.info(f"Compressing {out_vcf}")
    shell("pbgzip -c {out_vcf} > {snakemake.output['vcf']} 2> {log}")
    logging.info(f"Indexing {snakemake.output['call']}")
    shell("tabix -p vcf {snakemake.output['vcf']} >> {log} 2>&1")
    logging.info(f"Removing temporary file {out_vcf}")
    shell("rm --verbose {out_vcf} >> {log} 2>&1")