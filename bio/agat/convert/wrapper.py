"""Snakemake wrapper for agat conversion scripts"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2023, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"


from snakemake.shell import shell
log = snakemake.log_fmt_shell(stdout=True, stderr=True)
extra = snakemake.params.get("extra", "")

input_data = ""
output_data = ""
agat_script = snakemake.params.get("script")
outfile = snakemake.output[0]

# Configuration parameter is commont between all
# conversion sur-scripts
config = snakemake.input.get("config")
if config:
    input_data += f" --config {config} "


# Gathering input files
bed = snakemake.input.get("bed")
embl = snakemake.input.get("embl")
genscan = snakemake.input.get("genscan")
mfannot = snakemake.input.get("mfannot")
minimap2 = snakemake.input.get("minimap2")
gxf = snakemake.input.get("gxf")

# Only one script takes bed as input
if bed:
    input_data += f" --bed {bed} "
    output_data = f" --gff {outfile} "
    agat_script = agat_script or "agat_convert_bed2gff.pl "


# Only one script takes embl as input
elif embl:
    input_data += f" --embl {embl} "
    output_data = f" --gff {outfile} "
    agat_script = agat_script or "agat_convert_embl2gff.pl "


# Only one script takes genscan as input
elif genscan:
    input_data += f" --genscan {genscan} "
    output_data = f" --gff {outfile} "
    agat_script = agat_script or "agat_convert_genscan2gff.pl "


# Only one script takes mfannot as input
elif mfannot:
    input_data += f" --mfannot {mfannot} "
    output_data = f" --gff {outfile} "
    agat_script = agat_script or "agat_convert_mfannot2gff.pl "


# Only one script takes sam/bam/map as input
elif minimap2:
    if str(minimap2).lower().endswith(".sam"):
        input_data += f" --sam {minimap2} "
    elif str(minimap2).lower().endswith(".bam"):
        input_data += f" --bam {minimap2} "
    else:
        input_data += f" --input {minimap2} "
    output_data = f" --output {outfile} "
    agat_script = agat_script or "agat_convert_minimap2_bam2gff.pl "

# We need to know output file format and guess
# agat script to deduce gtf/gff/gxf input argument

# Only one script takes bed as output
elif str(outfile).lower().endswith(".bed"):
    input_data += f" --gff {gxf} "
    output_data = f" --outfile {outfile} "
    agat_script = agat_script or "agat_convert_sp_gff2bed.pl "


# Only one script takes gtf as output
elif str(outfile).lower().endswith("gtf"):
    # Possible ambiguity with agat_convert_sp_gxf2gxf.pl
    # Which also handles these file formats. However, special
    # gtf formats are handled by agat_convert_sp_gff2gtf.pl
    # and not by agat_convert_sp_gxf2gxf.pl
    # If user still wants to use agat_convert_sp_gxf2gxf.pl,
    # then they'll have to provide if in parameters.
    input_data += f" --gff {gxf} "
    output_data = f" --gtf {outfile} "
    agat_script = agat_script or "agat_convert_sp_gff2gtf.pl "


# Only one script takes tsv as output
elif str(outfile).lower().endswith("tsv"):
    input_data += f" --gff {gxf} "
    output_data = f" --outfile {outfile} "
    agat_script = agat_script or "agat_convert_sp_gff2tsv.pl "


# Only one script takes zff as output
elif str(outfile).lower().endswith("zff"):
    input_data += f" --gff {gxf} "
    output_data = f" --outfile {outfile} "
    agat_script = agat_script or "agat_convert_sp_gff2zff.pl "


# Only one script takes tsv as output
elif str(outfile).lower().endswith("gff"):
    input_data += f" --gff {gxf} "
    output_data = f" --outfile {outfile} "
    agat_script = agat_script or "agat_convert_sp_gxf2gxf.pl "


shell("{agat_script} {extra} {input_data} {output_data} {log}")