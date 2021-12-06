"""Snakemake wrapper for GATK4 Mutect2"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2019, Dayris Thibault"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

from snakemake.shell import shell
from snakemake.utils import makedirs
from snakemake_wrapper_utils.java import get_java_opts

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

bam_output = "--bam-output"
if snakemake.output.get("bam", None) is not None:
    bam_output = bam_output + " " + snakemake.output.bam
else:
    bam_output = ""


germline_resource = ""
if "germline" in snakemake.input.keys():
    germline_resource = "--germline-resource {}".format(
        snakemake.input["germline"]
    )


intervals = ""
if "intervals" in snakemake.input.keys():
    intervals = "--intervals {}".format(snakemake.input["intervals"])


f1r2 = ""
if "f1r2" in snakemake.output.keys():
    f1r2 = "--f1r2-tar-gz {}".format(snakemake.output["f1r2"])

tumor = ""
if "tumor" in snakemake.input.keys():
    tumor = "--input {}".format(snakemake.input["tumor"])

pon = ""
if "pon" in snakemake.input.keys():
    pon = "--panel-of-normals {}".format(snakemake.input["pon"])

extra = snakemake.params.get("extra", "")
java_opts = get_java_opts(snakemake)

java_opts += f" -XX:+UseParallelGC -XX:ParallelGCThreads={snakemake.threads}"

shell(
    "OMP_NUM_THREADS={snakemake.threads} && export OMP_NUM_THREADS && "
    "gatk --java-options '{java_opts}' Mutect2 "  # Tool and its subprocess
    " {tumor} "  # Path to tumor input file
    "--input {snakemake.input.map} "  # Path to input mapping file
    "{bam_output} "  # Path to output bam file, optional
    "{f1r2} "  # Path to output f1r2 count file
    "{germline_resource} "  # Path to optional germline resource VCF
    "{intervals} "  # Path to optional intervals
    "--output {snakemake.output.vcf} "  # Path to output vcf file
    "--reference {snakemake.input.fasta} "  # Path to reference fasta file
    "--native-pair-hmm-threads {snakemake.threads} "  # Maximum number of threads
    "{extra} "  # Extra parameters
    "{log}"  # Logging behaviour
)
