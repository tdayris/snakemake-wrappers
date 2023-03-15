#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np

germline_text = "Germline DNA"
ctc_text = "CTC"
wbc_text = "WBC"


def get_time_point(uid: str, tos: str) -> str:
    v = uid.strip().split("_")
    if tos.strip().lower() == germline_text.lower():
        return np.NaN
    if len(v) >= 2:
        if v[2].startswith("V") and len(v[2]) == 2:
            return float(v[2][1:])
        else:
            return np.NaN
    return np.NaN


def get_tp_treatment(tp: str, ttm: bool = False) -> str:
    tp = tp.strip().split(" - ")
    if ttm:
        if tp != ["Baseline"]:
            return tp[1]
        else:
            return "None"

    if tp[0] == "Baseline":
        return np.NaN
    return float(tp[0][1:])


def get_nb(uid: str, tos: str, target: bool = False) -> str:
    if tos.strip().lower() != germline_text.lower():
        if target:
            try:
                return float(uid.strip().split(" ")[-1])
            except ValueError:
                return np.NaN
        uid = "_".join(uid.strip().split(":")[0].split("_")[:-1])
        try:
            return float(uid.split("_")[-1])
        except ValueError:
            return np.NaN
    return np.NaN


def get_key(p: str, tp: float, tos: str, nb: float) -> str:
    return str(p) + "_V" + str(tp) + "_" + str(tos) + " " + str(nb)


def read_big_table(path: str) -> pd.DataFrame:
    print(f"Loading BigTableAnnoteted.tsv{path}")
    df = pd.read_csv(
        path, 
        sep=",", 
        header=None,
        low_memory=False
    )

    df.columns = [
        "CHROM",
        "POS",
        "REF",
        "ALT",
        "Method",
        "Gene_Name",
        "VARIANT_CLASS",
        "Consequence",
        "IMPACT",
        "Protein_position",
        "Amino_acids",
        "Hotspot",
        "Patient",
        "SampleName",
        "normal_ref_count",
        "tumor_ref_count",
        "normal_alt_count",
        "tumor_alt_count",
        "normal_depth",
        "tumor_depth",
        "normal_alt_freq",
        "tumor_alt_freq",
        "varID",
        "STRAND",
        "FILTER",
        "Feature_RefSeq",
        "CANONICAL",
        "Canonical_NM",
        "Canonical_Protein_position",
        "Canonical_Amino_acids",
        "Canonical_IMPACT",
        "Canonical_Consequence",
        "CODING",
        "inExon",
        "BIOTYPE",
        "Position_Type",
        "SIFT",
        "PolyPhen",
        "dbSNP_COSMIC",
        "uniqueID",
        "Sample_Type",
        "Tool",
        "Condition"
    ]

    chroms = list(range(1, 23)) + ["MT", "X", "Y"] + list(map(str, range(1, 23)))
    df = df[df["CHROM"].isin(chroms)]

    df["Type of sample"] = df["Condition"]

    print("Adding time points")
    df["Time-Point"] = [
        get_time_point(uid, tos) 
        for uid, tos in zip(df["uniqueID"], df["Type of sample"])
    ]

    print("Adding CTC number")
    df["CTC_nb"] = [
        get_nb(uid, tos) 
        for uid, tos in zip(df["uniqueID"], df["Type of sample"])
    ]

    df["KEY"] = [
        get_key(p, tp, tos, nb)
        for p, tp, tos, nb in zip(
            df["Patient"], df["Time-Point"], df["Type of sample"], df["CTC_nb"]
        )
    ]

    return df


def read_labels(path: str) -> pd.DataFrame:
    print(f"Loading {path}")
    df = pd.read_csv(path, sep="\t", header=0)

    print("Adding time points")
    df["Treatment"] = [
        get_tp_treatment(tp, True) for tp in df["Time-Point"]
    ]
    df["Time-Point"] = [
        get_tp_treatment(tp, False) for tp in df["Time-Point"]
    ]

    print("Adding CTC number")
    df["CTC_nb"] = [
        get_nb(uid, tos, True) 
        for uid, tos in zip(df["Targeted NGS SampleName"], df["Type of sample"])
    ]

    df["KEY"] = [
        get_key(p, tp, tos, nb)
        for p, tp, tos, nb in zip(
            df["Patient"], df["Time-Point"], df["Type of sample"], df["CTC_nb"]
        )
    ]

    df = df[[
        "Treatment",
        "FullComment",
        "Targeted NGS SampleName",
        "Number of cells",
        "Normal / Tumor",
        "WGA QC (/4)",
        "KEY"
    ]]

    return df
    

df = read_big_table(path=snakemake.input["bigtable"])
labels = read_labels(path=snakemake.input["label"]) # "EGFR_Targeted NGS_Sample sheet_21-07-2022-1.csv")

print(df.shape)
print(labels.shape)

merged = pd.merge(
    left=df,
    right=labels,
    how="left",
    on="KEY"
)

print(merged.shape)
merged.to_csv(snakemake.output["bigtable"], sep="\t", header=True, index=False)