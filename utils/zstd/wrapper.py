# coding: utf-8

from os.path import isdir, exists
from snakemake.shell import shell
from snakemake_wrapper_utils.snakemake import get_format, is_arg

log = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)
extra = snakemake.params.get("extra", "")
extra += " --force"

# Avoid automatic removal of input file which would lead to errors
# in snakemake IO controls:
if is_arg("--rm", extra):
    raise ValueError(
        "Input file should not be removed by `zstd` itself. "
        "Mark data as temporary with `temp` within Snakemake"
    )
    extra += " --keep"

outfiles = [str(o) for o in snakemake.output]
if len(outfiles) == 1:
    extra += f" -o {snakemake.output}"


# Guess compression format or decompression behavior
fmt = set([get_format(o, ignore_compression=False) for o in outfiles])
print(fmt)
if all(f == "gzip" for f in fmt):
    extra += " --format='gzip'"
elif all(f == "zst" for f in fmt):
    extra += " --format='zstd'"
elif all(f not in ("gzip", "zst") for f in fmt):
    extra += " --decompress"
else:
    raise ValueError(
        "Output files have divergent compression formats or "
        "mixed compressed/uncompressed files."
    )


shell("zstd -T{snakemake.threads} {extra} {snakemake.input} {log}")
