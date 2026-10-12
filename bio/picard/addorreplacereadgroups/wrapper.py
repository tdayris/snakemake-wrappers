__author__ = "Johannes Köster"
__copyright__ = "Copyright 2016, Johannes Köster"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"


import tempfile
from shlex import quote
from snakemake.shell import shell
from snakemake_wrapper_utils.java import get_java_opts

log = snakemake.log_fmt_shell()

extra = snakemake.params.get("extra", "")
java_opts = get_java_opts(snakemake)

bai = snakemake.output.get("bai")
if bai:
    extra += " --CREATE_INDEX"

md5 = snakemake.output.get("md5")
if md5:
    extra += " --CREATE_MD5_FILE"

ref = snakemake.input.get("ref")
if ref:
    extra += f" --REFERENCE_SEQUENCE {quote(ref)}"

with tempfile.TemporaryDirectory() as tmpdir:
    shell(
        "picard AddOrReplaceReadGroups"
        " {java_opts} {extra}"
        " --INPUT {snakemake.input.bam}"
        " --TMP_DIR {tmpdir}"
        " --OUTPUT {snakemake.output.bap}"
        " {log}"
    )
