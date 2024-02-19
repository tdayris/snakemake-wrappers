# Run CellRanger
"""
008.cellranger_arc_RNA_ATAC
from
-> 007.csv_cellranger_arc_RNA_ATAC
by
-> 003.multiqc
-> 003.irods_complient
"""

wildcard_constraints:
    CR_sample='|'.join([x for x in dic_DATA['singleCell_ATACseq_coupledScRNA'].keys()]),

#from 007.cellranger_csv_RNA_ATAC.smk
#def CR_fq_input_arc(wildcards):
#    return sum([ dic_DATA['singleCell_ATACseq_coupledScRNA'][wildcards.CR_sample][lib_type]['Fastq_Files'] for lib_type in list(dic_DATA['singleCell_ATACseq_coupledScRNA'][wildcards.CR_sample].keys()) ], [])

def ref_species_arc (wildcards):
    if dic_DATA['singleCell_ATACseq_coupledScRNA'][wildcards.CR_sample]["GE"]['Species'] == "human":
        return str('/mnt/beegfs/database/bioinfo/cellranger-arc/2.0.1/refdata-cellranger-arc-GRCh38-2020-A-2.0.0')
    elif dic_DATA['singleCell_ATACseq_coupledScRNA'][wildcards.CR_sample]["GE"]['Species'] == "mouse":
        return str('/mnt/beegfs/database/bioinfo/cellranger-arc/2.0.1/refdata-cellranger-arc-mm10-2020-A-2.0.0')
    else:
        sys.exit("Species not recognized for Gene Expression coupled to Chromatin analysis!")


rule cellranger_arc:
    input:
        CR_input=CR_fq_input_arc,
        csv_config="cellranger/CR_params_arc_{CR_sample}.csv"
    output:
        #temp("cellranger/web_summary/{CR_sample}_web_summary_mqc.html")
        temp("cellranger/csv_summary/{CR_sample}_metrics_summary.csv")
    message:
        "Computing CellRanger of {wildcards.CR_sample}"
    params:
        csv_config="CR_params_arc_{CR_sample}.csv",
        reference=ref_species_arc,
        sing_arg=workflow.singularity_args
    threads: 6
    resources:
        #mem_mb=lambda wildcard, attempt: min(attempt * 25000, 50000),
        mem_mb=64000,
        time_min=1430,
        #time_min=lambda wildcard, attempt: 2880,
        tmpdir="tmp"
    log:
        "logs/cellranger/{CR_sample}_cellranger_arc.log"
    shell:
        "module load singularity/3.6.3 && singularity exec --no-home {params.sing_arg} -B {PIPELINE_FOLDER}/scripts/:{PIPELINE_FOLDER}/scripts/ {PIPELINE_FOLDER}/envs/cellranger_arc_v2.0.2.simg {PIPELINE_FOLDER}/scripts/script_cellranger_arc_RNA_ATAC.sh {threads} {resources.mem_mb} {wildcards.CR_sample} {params.csv_config} {params.reference} $(pwd)"
