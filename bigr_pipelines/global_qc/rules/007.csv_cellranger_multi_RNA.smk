# Make csv config file for CellRanger
"""
007.cellranger_csv_RNA
from
-> Entry job
by
-> 008.cellranger_multi_RNA
"""

wildcard_constraints:
    CR_sample='|'.join([x for x in dic_DATA['SingleCell_RNAseq'].keys()]),


def CR_fq_input_RNA(wildcards):
    return sum([ dic_DATA['SingleCell_RNAseq'][wildcards.CR_sample][lib_type]['Fastq_Files'] for lib_type in list(dic_DATA['SingleCell_RNAseq'][wildcards.CR_sample].keys()) ], [])

def ref_species_ge (wildcards):
    if dic_DATA['SingleCell_RNAseq'][wildcards.CR_sample]["GE"]['Species'] == "human":
        return str('/mnt/beegfs/database/bioinfo/cellranger/6.0.2/refdata-gex-GRCh38-2020-A')
    elif dic_DATA['SingleCell_RNAseq'][wildcards.CR_sample]["GE"]['Species'] == "mouse":
        return str('/mnt/beegfs/database/bioinfo/cellranger/6.0.2/refdata-gex-mm10-2020-A')
    else:
        sys.exit("Species not recognized for Gene Expression analysis!")

def ref_species_vdj (wildcards):
    if ("TCR" in dic_DATA['SingleCell_RNAseq'][wildcards.CR_sample]):
        if (dic_DATA['SingleCell_RNAseq'][wildcards.CR_sample]["TCR"]['Species'] == "human"):
            return '/mnt/beegfs/database/bioinfo/cellranger/6.0.2/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0'
        elif (dic_DATA['SingleCell_RNAseq'][wildcards.CR_sample]["TCR"]['Species'] == "mouse"):
            return '/mnt/beegfs/database/bioinfo/cellranger/6.0.2/refdata-cellranger-vdj-GRCm38-alts-ensembl-5.0.0'
        else :
            sys.stderr.write("Species not recognized!")
    elif ("BCR" in dic_DATA['SingleCell_RNAseq'][wildcards.CR_sample]):
        if (dic_DATA['SingleCell_RNAseq'][wildcards.CR_sample]["BCR"]['Species'] == "human"):
            return '/mnt/beegfs/database/bioinfo/cellranger/6.0.2/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0'
        elif (dic_DATA['SingleCell_RNAseq'][wildcards.CR_sample]["BCR"]['Species'] == "mouse"):
            return '/mnt/beegfs/database/bioinfo/cellranger/6.0.2/refdata-cellranger-vdj-GRCm38-alts-ensembl-5.0.0'
        else :
            sys.exit("Species not recognized for VDJ analysis!")
    else: return ''
    
rule make_csv_for_cellranger_multi:
    input:
        CR_fq_input_RNA
    output:
        csv_config=temp("cellranger/CR_params_multi_{CR_sample}.csv")
    message:
        "Making csv parameters file of {wildcards.CR_sample}, for CellRanger analysis."
    params:
        ref_ge=ref_species_ge,
        ref_vdj=ref_species_vdj
    threads: 1
    resources:
        mem_gb=1,
        time_min=lambda wildcard, attempt: attempt * 5,
        tmpdir="tmp"
    log:
        "logs/cellranger/{CR_sample}_make_csv_for_cellranger_multi.log"
    run:
        ### Make params script
        f = open(output.csv_config, "w")
        
        #GE
        f.write("[gene-expression]\n")
        f.write("reference," + params.ref_ge + "\n")
        f.write("expect-cells,10000\n")
        f.write("no-secondary,false\n")
        f.write("no-bam,true\n")

        #VDJ
        if ("TCR" in dic_DATA['SingleCell_RNAseq'][wildcards.CR_sample] or "BCR" in dic_DATA['SingleCell_RNAseq'][wildcards.CR_sample]):
            f.write("[vdj]\n")
            f.write("reference," + params.ref_vdj + "\n")

        # Add libraries section
        #header
        f.write("[libraries]\n")
        f.write("fastq_id,fastqs,feature_types\n")
        #ge
        if "GE" in list(dic_DATA['SingleCell_RNAseq'][wildcards.CR_sample].keys()):
            for CR_sample_library in dic_DATA['SingleCell_RNAseq'][wildcards.CR_sample]['GE']['Library_Name']:
                #serach fastq files path
                fq_path = ''.join(frozenset([os.path.dirname(x) for x in dic_DATA['SingleCell_RNAseq'][wildcards.CR_sample]['GE']['Fastq_Files'] if CR_sample_library in x]))
                f.write(CR_sample_library + "," + str(fq_path) + ",gene expression\n")
        #tcr
        if "TCR" in list(dic_DATA['SingleCell_RNAseq'][wildcards.CR_sample].keys()):
            for CR_sample_library in dic_DATA['SingleCell_RNAseq'][wildcards.CR_sample]['TCR']['Library_Name']:
                #serach fastq files path
                fq_path = ''.join(frozenset([os.path.dirname(x) for x in dic_DATA['SingleCell_RNAseq'][wildcards.CR_sample]['TCR']['Fastq_Files'] if CR_sample_library in x]))
                f.write(CR_sample_library + "," + str(fq_path) + ",vdj-t\n")
        #bcr
        if "BCR" in list(dic_DATA['SingleCell_RNAseq'][wildcards.CR_sample].keys()):
            for CR_sample_library in dic_DATA['SingleCell_RNAseq'][wildcards.CR_sample]['BCR']['Library_Name']:
                #serach fastq files path
                fq_path = ''.join(frozenset([os.path.dirname(x) for x in dic_DATA['SingleCell_RNAseq'][wildcards.CR_sample]['BCR']['Fastq_Files'] if CR_sample_library in x]))
                f.write(CR_sample_library + "," + str(fq_path) + ",vdj-b\n")
        
        f.close()