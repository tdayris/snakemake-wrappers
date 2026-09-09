# coding: utf-8

"""Snakemake wrappers used process data with librarian locally"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2025, Thibault Dayris"
__license__ = "MIT"

from os import chdir, getcwd
from pathlib import Path
from snakemake.shell import shell
from snakemake_wrapper_utils.snakemake import move_files
from shlex import quote
from tempfile import TemporaryDirectory

extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

with TemporaryDirectory() as tempdir:
    librarian_directory = str(Path(snakemake.input.librarian_directory).resolve())
    launcher_script = f"{librarian_directory}/librarian_v1.1/librarian"
    reads = str(Path(snakemake.input.reads).resolve())

    expected_outputs = {
        "composition_map_svg": f"{tempdir}/librarian_compositions_map.svg",
        "composition_map_png": f"{tempdir}/librarian_compositions_map.png",
        "probability_maps_svg": f"{tempdir}/librarian_probability_maps.svg",
        "probability_maps_png": f"{tempdir}/librarian_probability_maps.png",
        "prediction_plot_svg": f"{tempdir}/librarian_prediction_plot.svg",
        "prediction_plot_png": f"{tempdir}/librarian_prediction_plot.png",
        "heatmap_txt": f"{tempdir}/librarian_librarian_heatmap.txt",
        "report": f"{tempdir}/librarian_Librarian_analysis.html",
    }

    move_cmds = " && ".join(move_files(snakemake, expected_outputs))

    shell(
        "( cd {tempdir} && "
        "  {launcher_script} --local {extra} {reads:q} && "
        "  ls {tempdir} && tree {tempdir} && "
        "  cd - && "
        "  {move_cmds} ) {log} "
    )
