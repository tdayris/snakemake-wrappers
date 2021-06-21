This table contains the results of the differential gene expression made with DESeq2, for `{{snakemake.wildcards.comparison}}`. It has been manually annotated, filtered and sorted. Only significative changes are kept.

+----------------+-----------------------------+
| Column name    | Content                     |
+================+=============================+
| Stat change    | The fold change value       |
+----------------+-----------------------------+
| P_Value        | The adjusted pvalue         |
+----------------+-----------------------------+
| GeneIdentifier | The unique name of the gene |
+----------------+-----------------------------+
| Gene_Name      | The human readable name     |
+----------------+-----------------------------+
| Chromosome     | The chromosome name         |
+----------------+-----------------------------+
| Start          | The gene starting position  |
+----------------+-----------------------------+
| Stop           | The gene end position       |
+----------------+-----------------------------+
| Stran          | The gene strand             |
+----------------+-----------------------------+
| cluster        | Up/Down regulation status   |
+----------------+-----------------------------+
| significance   | Significativity             |
+----------------+-----------------------------+

If a gene is labelled "Non_Significative" it does not mean a gene expression does not change. It means that the observed change has a too high risk of being driven by something else than the considered factor.

If a gene is missing, then it has been filtered out based on either its P-Value or its Fold Change. It is also possible that no expression has been measured on this gene and, consequently, it has been filtered out. Finally, it is possible that you use a common homonym for this gene, search its unique gene identifier.

This is a HTML file. This means you can open this file with your favorite browser, like Firefox, Chromium, Brave, etc. This HTML file is a local file, you do not need any Internet connection to open this file. In the same way, it will be "available" as long as you have this file on your computer. Do not erase the files and directories aside. These data are required for quick and efficient data display.
