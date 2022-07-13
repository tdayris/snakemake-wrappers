#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Format MANE summary to satisfy VCFTools requirements"""

import gzip
import logging

logging.basicConfig(
    filename=snakemake.log[0],
    filemode="w",
    level=logging.DEBUG
)

gunzipped_output = str(snakemake.output["mane"])[:-len(".gz")]

def annot(NCBI_GeneID: str,
          Ensembl_Gene: str,
          HGNC_ID: str,
          symbol: str,
          name: str,
          RefSeq_nuc: str,
          RefSeq_prot: str,
          Ensembl_nuc: str,
          Ensembl_prot: str,
          MANE_status: str,
          GRCh38_chr: str,
          chr_start: int,
          chr_end: int,
          chr_strand: str) -> str:
    return (
        f"MANE_NCBI_GeneID={NCBI_GeneID};"
        f"MANE_Ensembl_Gene={Ensembl_Gene};"
        f"MANE_HGNC_ID={HGNC_ID};"
        f"MANE_symbol={symbol};"
        f"MANE_name={name.replace(' ', '_').replace('-', '_').replace('/', '_')};"
        f"MANE_RefSeq_nuc={RefSeq_nuc};"
        f"MANE_RefSeq_prot={RefSeq_prot};"
        f"MANE_Ensembl_nuc={Ensembl_nuc};"
        f"MANE_Ensembl_prot={Ensembl_prot};"
        f"MANE_MANE_status={MANE_status};"
        f"MANE_GRCh38_chr={GRCh38_chr};"
        f"MANE_chr_start={chr_start};"
        f"MANE_chr_end={chr_end}"
    )


with gzip.open(snakemake.input["mane"], "r") as manin, gzip.open(snakemake.output["mane"], "wb") as manout:
    for line in manin:
        line = line.decode("utf-8")

        if line.startswith("#"):
            continue

        chomp = line[:-1].split("\t")

        try:
            new_line = "\t".join([
                chomp[10],
                chomp[11],
                chomp[12],
                annot(*chomp)
            ])
        except IndexError:
            logging.debug(chomp)
            raise

        new_line = bytes(new_line + "\n", "utf-8")
        manout.write(new_line)
