__author__ = "Christopher Schröder"
__copyright__ = "Copyright 2020, Christopher Schröder"
__email__ = "christopher.schroeder@tu-dortmund.de"
__license__ = "MIT"


from os import path

from snakemake.shell import shell


# Extract arguments.
extra = snakemake.params.get("extra", "")
sort_extra = snakemake.params.get("sort_extra", "")
samblaster_extra = snakemake.params.get("samblaster_extra", "")
sambamba_view_extra = snakemake.params.get("sambamba_view_extra", "")

index = snakemake.input.get("index", "")
if isinstance(index, str):
    index = path.splitext(snakemake.input.idx)[0]
else:
    index = path.splitext(snakemake.input.idx[0])[0]

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# Check inputs/arguments.
if not isinstance(snakemake.input.reads, str) and len(snakemake.input.reads) not in {
    1,
    2,
}:
    raise ValueError("input must have 1 (single-end) or " "2 (paired-end) elements")


threads = snakemake.threads
if threads < 4:
    raise ValueError("At least 4 threads are required for this wrapper")
if (threads-1)%3 == 0:
    bwa_threads = sambamba_view_threads = sambamba_sort_threads = int((threads-1)/3)
if (threads-1)%3 == 1:
    bwa_threads = int((threads-1)/3)+1
    sambamba_view_threads = sambamba_sort_threads = int((threads-1)/3)
if (threads-1)%3 == 2:
    bwa_threads = sambamba_sort_threads = int((threads-1)/3)+1
    sambamba_view_threads = int((threads-1)/3)

shell(
    "(bwa mem"
    " -t {bwa_threads}"
    " {extra}"
    " {index}"
    " {snakemake.input.reads}"
    " | samblaster"
    " {samblaster_extra}"
    " | sambamba view -S -f bam /dev/stdin"
    " -t {sambamba_view_threads}"
    " {sambamba_view_extra}"
    " | sambamba sort /dev/stdin"
    " -t {sambamba_sort_threads}"
    " -o {snakemake.output.bam}"
    " {sort_extra}"
    ") {log}"
)
