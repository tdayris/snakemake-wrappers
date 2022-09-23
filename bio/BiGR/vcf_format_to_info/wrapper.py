"""Snakemake wrapper which copies FORMAT to INFO"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2021, Dayris Thibault"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

from __future__ import annotations
import datetime
import gzip
import logging
import os.path

from snakemake.utils import makedirs
from typing import List, Optional, Union

logging.basicConfig(
    filename=snakemake.log[0],
    filemode="w",
    level=logging.DEBUG
)

makedirs(os.path.dirname(snakemake.output["call"]))


def open_function(file: str):
    """Return the correct opening function"""
    if file.endswith(".gz"):
        return gzip.open(file, "rb")
    return open(file, "r")


def get_supplementary_headers(prefixes: List[str]) -> str:
    prefixes = [p.replace("-", "_") for p in prefixes]
    headers = [
        """##FILTER=<ID=IsGermline,Number=.,Type=String,Description="Variant exists in Normal">\n""",
        """##FILTER=<ID=IsSomatic,Number=.,Type=String,Description="Variant does not exists in Normal, but exists in Tumor">\n""",
        """##INFO=<ID=Match_Norm_Seq_Allele1,Number=.,Type=String,Description="Normal allele 1">\n""",
        """##INFO=<ID=Match_Norm_Seq_Allele2,Number=.,Type=String,Description="Normal allele 2">\n""",
        """##INFO=<ID=Start_Position,Number=.,Type=String,Description="Mutation start coordinate (1 based)">\n""",
        """##INFO=<ID=End_Position,Number=.,Type=String,Description="Mutation end coordinate (1 based)">\n""",
        """##INFO=<ID=Tumor_Seq_Allele1,Number=.,Type=String,Description="Tumor allele 1">\n""",
        """##INFO=<ID=Tumor_Seq_Allele2,Number=.,Type=String,Description="Tumor allele 2">\n""",
        """##INFO=<ID=t_depth,Number=.,Type=String,Description="Read depth across this locus in tumor">\n""",
        """##INFO=<ID=t_ref_count,Number=.,Type=String,Description="Read depth supporting the reference allele in tumor">\n""",
        """##INFO=<ID=t_alt_count,Number=.,Type=String,Description="Read depth supporting the variant allele in tumor">\n""",
        """##INFO=<ID=n_depth,Number=.,Type=String,Description="Read depth across this locus in normal">\n""",
        """##INFO=<ID=n_ref_count,Number=.,Type=String,Description="Read depth supporting the reference allele in normal">\n""",
        """##INFO=<ID=n_alt_count,Number=.,Type=String,Description="Read depth supporting the variant allele in normal">\n""",
        """##INFO=<ID=vcf_tumor_gt,Number=.,Type=String,Description="Tumor sample genotype column from VCF">\n""",
        """##INFO=<ID=vcf_normal_gt,Number=.,Type=String,Description="Normal sample genotype column from VCF">\n""",
        """##INFO=<ID=MutationNumberAt,Number=.,Type=String,Description="This is the nth mutation seen at this position">\n""",
        """##INFO=<ID=Matched_Norm_Sample_Barcode,Number=.,Type=String,Description="Normal sample name">\n""",
        """##INFO=<ID=Tumor_Sample_Barcode,Number=.,Type=String,Description="Tumor sample name">\n""",
    ]

    for prefix in prefixes:
        headers += [
            f"""##INFO=<ID={prefix}_MutationStatus,Number=.,Type=String,Description="{prefix} mutation status (het, hom, ...)">\n""",
            f"""##INFO=<ID={prefix}_Reference_Allele,Number=.,Type=String,Description="{prefix} reference allele">\n""",
            f"""##INFO=<ID={prefix}_Seq_Allele1,Number=.,Type=String,Description="{prefix} alternative allele 1">\n""",
            f"""##INFO=<ID={prefix}_Seq_Allele2,Number=.,Type=String,Description="{prefix} alternative allele 2">\n""",
            f"""##INFO=<ID={prefix}_DP,Number=.,Type=String,Description="{prefix} read depth">\n""",
            f"""##INFO=<ID={prefix}_AD_allele1,Number=.,Type=String,Description="{prefix} allele 1 depth">\n""",
            f"""##INFO=<ID={prefix}_AD_allele2,Number=.,Type=String,Description="{prefix} allele 2 depth">\n""",
            f"""##INFO=<ID={prefix}_AF,Number=.,Type=String,Description="{prefix} allele frequency">\n""",
            f"""##INFO=<ID={prefix}_Seq_Allele2,Number=.,Type=String,Description="{prefix} alternative allele 2">\n"""
        ]

    return headers


def annotate_mutect2(genotype: str,
                     ref: str,
                     alt: str,
                     prefix: str,
                     filter: str,
                     dp: str,
                     ad: str,
                     af: str,
                     normal: bool = False,
                     tumor: bool = False) -> List[Union[List[str], str]]:
    """Break down a genotype information, break down DP/AD/AF into info"""

    prefix = prefix.replace("-", "_")
    annotation = [
        f"{prefix}_Reference_Allele={ref}",
        f"{prefix}_DP={dp}",
        f"{prefix}_AF={af}"
    ]

    if normal is True:
        annotation.append(f"Matched_Norm_Sample_Barcode={prefix}")
        annotation.append(f"n_depth={dp}")
        annotation.append(f"vcf_normal_gt={genotype}")
    if tumor is True:
        annotation.append(f"Tumor_Sample_Barcode={prefix}")
        annotation.append(f"t_depth={dp}")
        annotation.append(f"vcf_tumor_gt={genotype}")

    if genotype in ["./.", "././.", "./././."]:
        return annotation

    genotype_list = genotype.split("/")
    if genotype_list[0] == "0":
        annotation.append(f"{prefix}_Seq_Allele1={ref}")
        annotation.append(f"{prefix}_AD_allele1={ad.split(',')[0]}")
        if normal is True:
            annotation.append(f"Match_Norm_Seq_Allele1={ref}")
            annotation.append(f"n_ref_count={ad.split(',')[0]}")
        if tumor is True:
            annotation.append(f"Tumor_Seq_Allele1={ref}")
            annotation.append(f"t_ref_count={ad.split(',')[0]}")
    else:
        annotation.append(f"{prefix}_Seq_Allele1={alt}")
        if normal is True:
            annotation.append(f"Match_Norm_Seq_Allele1={alt}")
        if tumor is True:
            annotation.append(f"Tumor_Seq_Allele1={alt}")

    for idx, value in enumerate(genotype_list[1:], start=2):
        try:
            annotation.append(f"{prefix}_AD_allele{idx}={ad.split(',')[idx-1]}")
            if value == "0":
                annotation.append(f"{prefix}_Seq_Allele{idx}={ref}")
                if normal is True:
                    annotation.append(f"Match_Norm_Seq_Allele{idx}={ref}")
                if tumor is True:
                    annotation.append(f"Tumor_Seq_Allele{idx}={ref}")
            else:
                    annotation.append(f"{prefix}_Seq_Allele{idx}={alt}")
                    annotation.append(f"t_alt_count={ad.split(',')[idx-1]}")
                    if normal is True:
                        annotation.append(f"Match_Norm_Seq_Allele{idx}={alt}")
                        annotation.append(f"n_alt_count={ad.split(',')[idx-1]}")
                    if tumor is True:
                        annotation.append(f"Tumor_Seq_Allele{idx}={alt}")
                        annotation.append(f"MutationNumberAt={idx-1}")
        except IndexError:
            logging.error(
                f"Genotype has a different length compared to corresponding allele depth. Skipping: {genotype}, {prefix}, {ad}"
            )

    if all(i == "0" for i in genotype_list):
        annotation.append(f"{prefix}_MutationStatus=HomozygousReference")
    elif all(i == "1" for i in genotype_list):
        annotation.append(f"{prefix}_MutationStatus=HomozygousMutant")
    else:
        annotation.append(f"{prefix}_MutationStatus=Heterozygous")
    logging.debug(f"{filter}, {str(annotation)}")
    return filter or "", ";".join(annotation)

colnames = None
samples = None
version = 2.0
name = "vcf_format_to_info"
url = f"github.com/tdayris/snakemake-wrappers/tree/Unofficial/bio/BiGR/{name}/wrapper.py"
headers = [f'##BiGRCommandLine=<ID={name},CommandLine="{url}",Version="{version}",Date="{datetime.date.today()}">\n']


if str(snakemake.output["call"]).endswith("vcf.gz"):
    out_vcf = snakemake.output["call"][:-3]
else:
    out_vcf = snakemake.output["call"]

with (open_function(snakemake.input["call"]) as instream,
      open(out_vcf, "w") as outstream):

    for line in instream:
        if isinstance(line, bytes):
            line = line.decode("utf-8")

        if line.startswith("##"):
            outstream.write(line)
            continue
        if line.startswith("#"):
            colnames = line[1:-1].split("\t")
            samples = colnames[9:]

            headers += get_supplementary_headers(samples)
            outstream.write("".join(headers))
            outstream.write(line)
        else:
            chomp = dict(zip(colnames, line[:-1].split("\t")))

            for idx, sample in enumerate(samples):
                format_sample = dict(zip(chomp["FORMAT"].split(":"), chomp[sample].split(":")))
                new_annotation = annotate_mutect2(
                    genotype = format_sample.get("GT", "./."),
                    ref = chomp["REF"],
                    alt = chomp["ALT"],
                    prefix = sample,
                    filter = chomp["FILTER"],
                    ad = format_sample.get("AD", ".,."),
                    af = format_sample.get("AF", "."),
                    dp = format_sample.get("DP", "."),
                    normal = snakemake.params.get("normal_sample", None) == sample,
                    tumor = snakemake.params.get("tumor_sample", None) == sample
                )
                if new_annotation[0] is not None:
                    if new_annotation != chomp["FILTER"]:
                        chomp["FILTER"] = new_annotation[0]

                if chomp["INFO"] in ["." or ""]:
                    chomp["INFO"] = new_annotation[1]
                else:
                    chomp["INFO"] += f";{new_annotation[1]}"


            outstream.write("\t".join([chomp[col] for col in colnames]) + "\n")


if str(snakemake.output["call"]).endswith("vcf.gz"):
    logging.info(f"Compressing {out_vcf}")
    shell("pbgzip -c {out_vcf} > {snakemake.output['vcf']} 2> {log}")
    logging.info(f"Indexing {snakemake.output['call']}")
    shell("tabix -p vcf {snakemake.output['vcf']} >> {log} 2>&1")
    logging.info(f"Removing temporary file {out_vcf}")
    shell("rm --verbose {out_vcf} >> {log} 2>&1")
