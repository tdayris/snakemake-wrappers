#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""This is the Snakemake Wrapper to extract mutations from mixcr"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2020, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"


def parse_d(mut: str):
    """
    Parse deletions
    """
    for key, val in enumerate(mut):
        if val.isdigit():
            return {
                "type": "Deletion",
                "Pos": mut[key:],
                "Ref": mut[1:key],
                "Alt": "-"
            }

def parse_i(mut: str):
    """
    Parse insertions
    """
    for key, val in enumerate(mut):
        if not val.isdigit() and val != "I":
            return {
                "type": "Insertion",
                "Pos": mut[1:key],
                "Ref": "-",
                "Alt": mut[key:]
            }


def parse_s(mut: str):
    """
    Parse substitutions
    """
    ref_limit = None
    alt_limit = None
    for key, val in enumerate(mut):
        if val.isdigit():
            ref_limit = key
            break
    for key, val in sorted(enumerate(mut), reverse=True):
        if val.isdigit():
            alt_limit = key
            break
    return {
        "type": "Substitution",
        "Pos": mut[ref_limit:alt_limit+1],
        "Ref": mut[1:ref_limit],
        "Alt": mut[alt_limit+1:]
    }


def parse_mutations(mutation_string: str):
    """
    Parse a mutation string into dict of dicts
    """
    separator = {"D": parse_d, "S": parse_s, "I": parse_i}

    if mutation_string == "":
        return [""]

    current_mut = [None, None]
    for val in mutation_string:
        if val in separator:
            if all(i is not None for i in current_mut):
                yield separator[current_mut[0]](current_mut[1])
                current_mut = [None, None]

            current_mut[0] = val
            current_mut[1] = val
        else:
            current_mut[1] += val
    yield separator[current_mut[0]](current_mut[1])


column_to_parse = snakemake.params.get("column", "allVAlignments")
gene_feature_name = snakemake.params.get("gene", "bestVHitScore")

with open(snakemake.input[0]) as infile, open(snakemake.output[0], "w") as outfile:
    lines = iter(infile)
    header = next(lines)

    positions = [None, None]
    for pos, val in enumerate(header[:-1].split("\t")):
        if val == gene_feature_name:
            positions[0] = pos - 1
        elif val == column_to_parse:
            positions[1] = pos

        if all(var is not None for var in positions):
            break
    else:
        raise ValueError("Could not find mutations, or feature names")

    line = next(lines)
    while line:
        chomp = line[:-1].split("\t")
        feature = chomp[positions[0]]
        mutations_string = chomp[positions[1]]
        mutation_list = [
            i.split("|")[-2] for i in mutations_string.split(",") if i != ""
        ]
        for mutations in mutation_list:
            for mutation in parse_mutations(mutations):
                if mutation is not None:
                    newline = "\t".join([
                        feature,
                        mutation["Pos"],
                        mutation["Ref"],
                        mutation["Alt"]
                    ])
                    print(newline)
                    outfile.write(newline + "\n")
        line = next(lines, None)