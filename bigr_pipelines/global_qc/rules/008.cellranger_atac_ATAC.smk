# Run CellRanger
"""
008.cellranger_atac_ATAC
from
-> Entry job
by
-> 003.multiqc
-> 003.irods_complient
"""

wildcard_constraints:
    CR_sample='|'.join([x for x in dic_DATA['singleCell_ATACseq'].keys()]),


def CR_library_names_ATAC(wildcards):
    return str(",".join(dic_DATA['singleCell_ATACseq'][wildcards.CR_sample]['ATAC']['Library_Name']))

def CR_fq_input_ATAC(wildcards):
    return dic_DATA['singleCell_ATACseq'][wildcards.CR_sample]['ATAC']['Fastq_Files']

def CR_fq_path_ATAC(wildcards):
    return str(','.join(frozenset([os.path.dirname(x) for x in dic_DATA['singleCell_ATACseq'][wildcards.CR_sample]['ATAC']['Fastq_Files']])))

def ref_species_atac (wildcards):
    if dic_DATA['singleCell_ATACseq'][wildcards.CR_sample]["ATAC"]['Species'] == "human":
        return str('/mnt/beegfs/database/bioinfo/cellranger-atac/2.1.0/refdata-cellranger-arc-GRCh38-2020-A-2.0.0')
    elif dic_DATA['singleCell_ATACseq'][wildcards.CR_sample]["ATAC"]['Species'] == "mouse":
        return str('/mnt/beegfs/database/bioinfo/cellranger-atac/2.1.0/refdata-cellranger-arc-mm10-2020-A-2.0.0')
    else:
        sys.exit("Species not recognized for Gene Expression coupled to Chromatin analysis!")


rule cellranger_atac:
    input:
        CR_input=CR_fq_input_ATAC
    output:
        #temp("cellranger/web_summary/{CR_sample}_web_summary_mqc.html")
        temp("cellranger/csv_summary/{CR_sample}_metrics_summary.csv")
    message:
        "Computing CellRanger of {wildcards.CR_sample}"
    params:
        reference=ref_species_atac,
        library_names=CR_library_names_ATAC,
        path_fastqs=CR_fq_path_ATAC,
        sing_arg=workflow.singularity_args
    threads: 6
    resources:
        #mem_mb=lambda wildcard, attempt: min(attempt * 25000, 50000),
        mem_mb=64000,
        time_min=1430,
        #time_min=lambda wildcard, attempt: 2880,
        tmpdir="tmp"
    log:
        "logs/cellranger/{CR_sample}_cellranger_atac.log"
    shell:
        "module load singularity/3.6.3 && singularity exec --no-home {params.sing_arg} -B {PIPELINE_FOLDER}/scripts/:{PIPELINE_FOLDER}/scripts/ {PIPELINE_FOLDER}/envs/cellranger_atac_v2.1.0.simg {PIPELINE_FOLDER}/scripts/script_cellranger_atac_ATAC.sh {threads} {resources.mem_mb} {wildcards.CR_sample} {params.library_names} {params.reference} {params.path_fastqs} $(pwd)"
