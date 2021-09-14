"""Snakemake wrapper to annotate variant occurence"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2021, Dayris Thibault"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"


import bgzip
import datetime
import logging

logging.basicConfig(
    filename=snakemake.log[0],
    filemode="w",
    level=logging.INFO
)

def annotate_vcf(inv: str, outv: str, occurence: dict[str, int]) -> None:
    """
    Use occurence to annotate VCF
    """
    logging.info(f"Working on {inv}")
    with open(inv, "r") as invcf, open(outv, 'w') as outvcf:
        for line in invcf:
            if line.startswith("##"):
                # Base header, nothing to do
                pass

            elif line.startswith("#"):
                # End of header, our new header must be added
                outvcf.write(header)

            else:
                chomp = line.split("\t")
                var_occ = occurence[":".join(chomp[0:5])]
                if chomp[7] == ".":
                    chomp[7] = f"VarOcc={var_occ}"
                else:
                    chomp[7] += f";VarOcc={var_occ}"

                line = "\t".join(chomp)

            outvcf.write(line)


def build_occurence_dict(occ_path: str) -> dict[str, int]:
    """
    From a text file, build the occurence dict
    """
    occurences = {}
    with open(occ_path) as occ_stream:
        for line in occ_stream:
            occ, chr, pos, id, ref, alt = line[:-1].split("\t")
            occurences[f"{chr}:{pos}:{id}:{ref}:{alt}"] = occ
    logging.debug("Variant occurences within cohort loaded")
    return occurences

version = 1.0
name = "VariantOccurence/Annotate"
url = f"github.com/tdayris/snakemake-wrappers/tree/Unofficial/bio/variantoccurence/{name}/wrapper.py"
header = f"""##BiGRCommandLine=<ID={name},CommandLine="{url}",Version={version},Date={datetime.date.today()}>\n"""
print("helloworld")
occ_dict = build_occurence_dict(snakemake.input["occurence"])
for incall, outcall in zip(snakemake.input["calls"], snakemake.output["calls"]):
    annotate_vcf(incall, outcall, occ_dict)
