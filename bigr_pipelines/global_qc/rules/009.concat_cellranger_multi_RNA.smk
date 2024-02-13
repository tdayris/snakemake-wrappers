# Run CellRanger
"""
009.concat_cellranger_multi_RNA
from
-> 008.cellranger_multi_RNA
by
-> 003.multiqc
-> 003.irods_complient
"""

wildcard_constraints:
    CR_sample='|'.join([x for x in dic_DATA['SingleCell_RNAseq'].keys()]),


def CR_csv_input_RNA(wildcards):
    if isScRNAData:
        all_samples = [x for x in dic_DATA['SingleCell_RNAseq'].keys()]
        return ["cellranger/csv_summary/" + s + "_metrics_summary.csv" for s in all_samples ]
    else:
        return []


rule concat_cellranger_multi:
    input:
        CR_input=CR_csv_input_RNA
    output:
        concat_multi=temp("cellranger/csv_concat/CellRanger_RNA_summary_mqc.csv")
    message:
        "Computing concatenation of CelleRanger multi metrics_summary.csv files"
    threads: 1
    resources:
        mem_mb=lambda wildcard, attempt: attempt * 1000,
        time_min=lambda wildcard, attempt: 60,
        tmpdir="tmp"
    log:
        "logs/cellranger/concat_cellranger_multi.log"
    run:
        #création d'un df vide
        final_df = pd.DataFrame()
        
        #pour chaque fichier
        for csv_files in input.CR_input:
            
            #lire le fichier
            tmp_df = pd.read_csv(csv_files)
            
            #ne sélectionner que les lignes dont la première colonne contient "Cells" (les autres métriques sont superflues et complexes à intégrer)
            tmp_df = tmp_df.loc[tmp_df['Category'] == 'Cells']
            
            #supprimer les colonnes "Category", "Grouped By" et "Group Name"
            tmp_df = tmp_df.drop(columns=['Category', 'Grouped By', 'Group Name'])
            
            #fusionner les colonnes "Library Type" et "Metric Name" (ex: "Cells: Confidently mapped reads in cells"), puis les supprimer, et ordonner les colonnes restantes
            tmp_df["Metric"] = tmp_df.apply(lambda x: ': '.join(x[['Library Type', 'Metric Name']]),axis=1)
            tmp_df = tmp_df.drop(columns=['Library Type', 'Metric Name'])
            tmp_df = tmp_df[["Metric", "Metric Value"]]
            
            #transposer le tableau et première ligne en tant que noms de colonnes.
            tmp_df = tmp_df.transpose()
            tmp_df.columns = tmp_df.iloc[0]
            tmp_df = tmp_df.drop(tmp_df.index[0])
            
            #ajouter le nom de l'échantillon en première colonne
            tmp_df.insert(0, "Sample ID", csv_files.replace('cellranger/csv_summary/', '').replace('_metrics_summary.csv', ''))

            #merger avec le tableau général
            final_df = pd.concat([final_df, tmp_df], ignore_index=True, sort=False)
        
        #order column
        final_df.insert(0, "Sample ID", final_df.pop("Sample ID"))
        final_df.insert(1, "Gene Expression: Cells", final_df.pop("Gene Expression: Cells"))
        final_df.insert(2, "Gene Expression: Median reads per cell", final_df.pop("Gene Expression: Median reads per cell"))
        
        #replace "," by " "
        final_df = final_df.replace(to_replace=r',', value='', regex=True)
        final_df.columns = [w.replace(',', ' and') for w in final_df.columns]
        final_df.columns = [w.replace('"', '') for w in final_df.columns]
        
        #sauvegarder le df en format csv
        final_df.to_csv(output.concat_multi, sep=',',header=True, index=False, na_rep='0', decimal='.')
        

        
        
        
        
        
        
        
        
        
        
        
        
        
        