# Run CellRanger
"""
009.concat_cellranger_arc_RNA_ATAC
from
-> 008.cellranger_arc_RNA_ATAC
by
-> 003.multiqc
-> 003.irods_complient
"""

wildcard_constraints:
    CR_sample='|'.join([x for x in dic_DATA['singleCell_ATACseq_coupledScRNA'].keys()]),


def CR_csv_input_RNA_ATAC(wildcards):
    if isScRNAATACData:
        all_samples = [x for x in dic_DATA['singleCell_ATACseq_coupledScRNA'].keys()]
        return ["cellranger/csv_summary/" + s + "_metrics_summary.csv" for s in all_samples ]
    else:
        return []


rule concat_cellranger_arc:
    input:
        CR_input=CR_csv_input_RNA_ATAC
    output:
        concat_arc=temp("cellranger/csv_concat/CellRanger_RNA_ATAC_summary_mqc.csv")
    message:
        "Computing concatenation of CellRanger arc summary.csv files"
    threads: 1
    resources:
        mem_mb=lambda wildcard, attempt: attempt * 1000,
        time_min=lambda wildcard, attempt: 60,
        tmpdir="tmp"
    log:
        "logs/cellranger/concat_cellranger_arc.log"
    run:
        #création d'un df vide
        final_df = pd.DataFrame()
        
        #pour chaque fichier
        for csv_files in input.CR_input:
            
            #lire le fichier
            tmp_df = pd.read_csv(csv_files)

            #merger avec le tableau général
            final_df = pd.concat([final_df, tmp_df], ignore_index=True, sort=False)
        
        #order column
        final_df.insert(0, "Sample ID", final_df.pop("Sample ID"))
        final_df.insert(1, "Estimated number of cells", final_df.pop("Estimated number of cells"))
        final_df.insert(2, "GEX Mean raw reads per cell", final_df.pop("GEX Mean raw reads per cell"))
        final_df.insert(3, "ATAC Mean raw read pairs per cell", final_df.pop("ATAC Mean raw read pairs per cell"))
        
        #replace "," by " "
        final_df = final_df.replace(to_replace=r',', value='', regex=True)
        
        #sauvegarder le df en format csv
        final_df.to_csv(output.concat_arc, sep=',',header=True, index=False, na_rep='NA', decimal='.')
        

        
        
        
        
        
        
        
        
        
        
        
        
        
        