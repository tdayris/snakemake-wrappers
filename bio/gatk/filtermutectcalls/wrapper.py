__author__ = "Patrik Smeds"
__copyright__ = "Copyright 2021, Patrik Smeds"
__email__ = "patrik.smeds@gmail.com"
__license__ = "MIT"


from snakemake.shell import shell
from snakemake_wrapper_utils.java import get_java_opts

extra = snakemake.params.get("extra", "")
java_opts = get_java_opts(snakemake)

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

bam = ""
if "bam" in snakemake.input.keys():
    bam = " --input {}".format(snakemake.input["bam"])


contamination = ""
if "contamination" in snakemake.input.keys():
    contamination = "--contamination-table {}".format(
        snakemake.input["contamination"]
    )

f1r2 = ""
if "f1r2" in snakemake.input.keys():
    f1r2 = "--orientation-bias-artifact-priors {}".format(
        snakemake.input["f1r2"]
    )

shell(
    "gatk --java-options '{java_opts}' FilterMutectCalls "
    "-R {snakemake.input.ref} -V {snakemake.input.vcf} "
    "{extra} {contamination} {f1r2} {bam} "
    "-O {snakemake.output.vcf} "
    "{log}"
)
