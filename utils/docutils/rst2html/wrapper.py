# coding: utf-8

"""Snakemake wrappers for rst2html"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2025, Thibault Dayris"
__license__ = "MIT"

from snakemake.shell import shell
from shlex import quote

extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

config = snakemake.input.get("config")
if config:
    extra += f" --config={quote(config)}"

warnings = snakemake.output.get("warnings")
if warnings:
    extra += f" --warnings={quote(warnings)}"

dependencies = snakemake.output.get("dependencies")
if dependencies:
    extra += f" --record-dependencies={quote(dependencies)}"

url_template = snakemake.input.get("url_template")
if url_template:
    extra += f" --pep-file-url-template={quote(url_template)}"

rfc_url = snakemake.input.get("rfc_url")
if rfc_url:
    extra += f" --rfc-base-url={quote(rfc_url)}"

html_template = snakemake.input.get("html_template")
if html_template:
    extra += f" --template={quote(html_template)}"

stylesheet = snakemake.input.get("stylesheet")
if stylesheet:
    if os.path.isdir(stylesheet):
        extra += f" --stylesheet-dirs={quote(stylesheet)}"
    else:
        extra += f" --stylesheet-path={quote(stylesheet)}"

shell(
    "rst2html {extra} --output-path='{snakemake.output.html}' "
    "{snakemake.input.rst} {log}"
)
