# Make csv config file for CellRanger
"""
007.cellranger_csv
from
-> Entry job
by
-> 008.cellranger_arc
"""

wildcard_constraints:
    CR_sample='|'.join([x for x in dic_DATA['singleCell_ATACseq_coupledScRNA'].keys()]),


def CR_fq_input_arc(wildcards):
    return sum([ dic_DATA['singleCell_ATACseq_coupledScRNA'][wildcards.CR_sample][lib_type]['Fastq_Files'] for lib_type in list(dic_DATA['singleCell_ATACseq_coupledScRNA'][wildcards.CR_sample].keys()) ], [])

rule make_csv_for_cellranger_arc:
    input:
        CR_fq_input_arc
    output:
        csv_config=temp("cellranger/CR_params_arc_{CR_sample}.csv")
    message:
        "Making csv parameters file of {wildcards.CR_sample}, for CellRanger analysis."
    threads: 1
    resources:
        mem_gb=1,
        time_min=lambda wildcard, attempt: attempt * 5,
        tmpdir="tmp"
    log:
        "logs/cellranger/{CR_sample}_make_csv_for_cellranger_arc.log"
    run:
        ### Make params script
        f = open(output.csv_config, "w")

        # Add libraries section
        #header
        f.write("fastqs,sample,library_type\n")

        #ge
        if "GE" in list(dic_DATA['singleCell_ATACseq_coupledScRNA'][wildcards.CR_sample].keys()):
            for CR_sample_library in dic_DATA['singleCell_ATACseq_coupledScRNA'][wildcards.CR_sample]['GE']['Library_Name']:
                #serach fastq files path
                fq_path = ''.join(frozenset([os.path.dirname(x) for x in dic_DATA['singleCell_ATACseq_coupledScRNA'][wildcards.CR_sample]['GE']['Fastq_Files'] if CR_sample_library in x]))
                f.write(str(fq_path) + "," + CR_sample_library + ",Gene Expression\n")
        #atac
        if "ATAC" in list(dic_DATA['singleCell_ATACseq_coupledScRNA'][wildcards.CR_sample].keys()):
            for CR_sample_library in dic_DATA['singleCell_ATACseq_coupledScRNA'][wildcards.CR_sample]['ATAC']['Library_Name']:
                #serach fastq files path
                fq_path = ''.join(frozenset([os.path.dirname(x) for x in dic_DATA['singleCell_ATACseq_coupledScRNA'][wildcards.CR_sample]['ATAC']['Fastq_Files'] if CR_sample_library in x]))
                f.write(str(fq_path) + "," + CR_sample_library + ",Chromatin Accessibility\n")

        f.close()