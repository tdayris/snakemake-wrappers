# Run CellRanger
"""
008.cellranger_multi
from
-> 007.cellranger_csv_RNA
by
-> 003.multiqc
-> 003.irods_complient
"""

wildcard_constraints:
    CR_sample='|'.join([x for x in dic_DATA['SingleCell_RNAseq'].keys()]),

#from 007.cellranger_csv_RNA.smk
#def CR_fq_input_RNA(wildcards):
#    return [ dic_DATA['SingleCell_RNAseq'][wildcards.CR_sample][lib_type]['Fastq_Files'] for lib_type in dic_DATA['SingleCell_RNAseq'][wildcards.CR_sample].keys() ]

rule cellranger_multi:
    input:
        CR_input=CR_fq_input_RNA,
        csv_config="cellranger/CR_params_multi_{CR_sample}.csv"
    output:
        #temp("cellranger/web_summary/{CR_sample}_web_summary_mqc.html")
        temp("cellranger/csv_summary/{CR_sample}_metrics_summary.csv")
    message:
        "Computing CellRanger of {wildcards.CR_sample}"
    params:
        csv_config="CR_params_multi_{CR_sample}.csv",
        sing_arg=workflow.singularity_args
    threads: 3
    resources:
        #mem_mb=lambda wildcard, attempt: min(attempt * 25000, 50000),
        mem_mb=25000,
        time_min=1430,
        #time_min=lambda wildcard, attempt: 2880,
        tmpdir="tmp"
    log:
        "logs/cellranger/{CR_sample}_cellranger_multi.log"
    shell:
        "module load singularity/3.6.3 && singularity exec --no-home {params.sing_arg} ../envs/cellranger_v7.2.0.simg ../scripts/script_cellranger_multi_RNA_TCR_BCR.sh {threads} {resources.mem_mb} {wildcards.CR_sample} {params.csv_config} $(pwd)"
