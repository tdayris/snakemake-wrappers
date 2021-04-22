__author__ = "Jan Forster"
__copyright__ = "Copyright 2020, Jan Forster"
__email__ = "j.forster@dkfz.de"
__license__ = "MIT"


from snakemake.shell import shell
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

if snakemake.threads != 2:
    raise ValueError("BCFTools reheader wrapper requires exactly 2 threads")

## Extract arguments
header = snakemake.input.get("header", None)
if header is not None:
    header_cmd = "--header " + header
else:
    header_cmd = ""

samples = snakemake.input.get("samples", None)
if samples is not None:
    samples_cmd = "--samples " + samples
else:
    samples_cmd = ""

fai = snakemake.input.get("fai", None)
if fai is not None:
    fai_cmd = "--fai " + fai
else:
    fai_cmd = ""


extra = snakemake.params.get("extra", "")
view_extra = snakemake.params.get("view_extra", "")

if str(snakemake.output).endswith(".gz"):
    view_extra += " --output-type z "
elif str(snakemake.output).endswith(".bcf"):
    view_extra += " --output-type b "
elif str(snakemake.output).endswith(".vcf"):
    view_extra += " --output-type v "

shell(
    "bcftools reheader "
    "{extra} "
    "{header_cmd} "
    "{samples_cmd} "
    "{fai_cmd} "
    "{snakemake.input.vcf} "
    "| bcftools view "
    "{view_extra} "
    "> {snakemake.output} "
    "{log}"
)
