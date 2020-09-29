__author__ = "Johannes Köster"
__copyright__ = "Copyright 2016, Johannes Köster"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"


from snakemake.shell import shell
log = snakemake.log_fmt_shell(stdout=True, stderr=True)


extra = snakemake.params[0]
output_file = snakemake.output[0]
if output_file.endswith(".bcf"):
    extra += " --output-type b "
elif output_file.endswith(".vcf.gz"):
    extra += " --output-type z "
elif output_file.endswith(".vcf"):
    extra += " --output-type v "
else:
    raise ValueError(
        "Output file extension should be one of: vcf, bcf, vcf.gz"
    )

shell(
    " bcftools view "
    " {extra} "
    " {snakemake.input[0]} "
    " -o {output_file} "
    " {log} "
)
