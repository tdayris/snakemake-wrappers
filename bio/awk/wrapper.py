#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Snakemake wrapper for Awk scripts"""


from snakemake.shell import shell
from typing import List, Optional

def join_awk_sections(instructions: Optional[List[str]] = None) -> Optional[str]:
    """Return the joint awk instruction, of None is no instruction are provided"""
    if (instructions == []) or (instructions is None):
        return None
    if isinstance(instructions, str):
        return instructions
    return ";".join(instructions)

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

begin = ""
if (begin_instructions := join_awk_sections(snakemake.params.get("begin", []))) is not None:
    begin = f'BEGIN{{{begin_instructions}}}'

end = ""
if (end_instructions := join_awk_sections(snakemake.params.get("end", []))) is not None:
    end = f'END{{{end_instructions}}}'

body = f'{{{join_awk_sections(snakemake.params.get("body", []))}}}'

shell(
    "awk '{begin} {body} {end}' {snakemake.input} > {snakemake.output} {log}"
)