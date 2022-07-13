__author__ = "Dayne Filer"
__copyright__ = "Copyright 2019, Dayne Filer"
__email__ = "dayne.filer@gmail.com"
__license__ = "MIT"


from snakemake.shell import shell
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

extra = snakemake.params.get("extra", "")
fasta = ""
if "fasta" in snakemake.input.keys():
    fasta = "--fasta-ref {}".format(snakemake.input["fasta"])

output_type = "--output-type "
if snakemake.output[0].endswith("bcf"):
    output_type += "b"
elif snakemake.output[0].endswith("vcf"):
    output_type += "v"
elif snakemake.output[0].endswith("vcf.gz"):
    output_type += "z"
else:
    output_type = "u"

region = ""
if "region" in snakemake.input.keys():
    region = "--regions {}".format(snakemake.input["region"])

shell(
    "bcftools norm {extra} {region} {fasta} "
    "{output_type} --threads {snakemake.threads} "
    "{snakemake.input.call} -o {snakemake.output[0]} {log}"
)
