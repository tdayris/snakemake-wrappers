This table contains the results of the differential gene expression made with DESeq2, for `{{snakemake.wildcards.comparison}}`. It has been manually annotated. Only differentially
expressed genes are in this file, and they are sorted by adjusted p-values.

This file contains the following information:

+-----------------+-----------------------------+
| Column name     | Content                     |
+=================+=============================+
| GeneIdentifier  | The unique name of the gene |
+-----------------+-----------------------------+
| Adjusted_PValue | The adjusted pvalue         |
+-----------------+-----------------------------+
| pvalue          | The raw pvalue              |
+-----------------+-----------------------------+
| stat_change     | The log2(fold change)       |
+----------------+-----------------------------+
| Gene_Name      | The human readable name     |
+----------------+-----------------------------+
| Chromosome     | The chromosome name         |
+----------------+-----------------------------+
| Start          | The gene starting position  |
+----------------+-----------------------------+
| Stop           | The gene end position       |
+----------------+-----------------------------+
| Strand         | The gene strand             |
+----------------+-----------------------------+
| cluster        | Up/Down regulation status   |
+----------------+-----------------------------+
| significance   | Significativity             |
+----------------+-----------------------------+
| ...            | Normalized counts           |
+----------------+-----------------------------+

If a gene is labelled "Non_Significative" it does not mean a gene expression does not change. It means that the observed change has a too high risk of being driven by something else than the considered factor.

If a gene is missing, then it is possible that no expression has been measured on this gene and, consequently, it has been filtered out. Finally, it is possible that you use a common homonym for this gene, search its unique gene identifier.

This is a TSV-formatted file. It can be opened in your favorite tabular file reader, like LibreOffice Calc, Excel, etc. Open your tabular file reader, then hit "open file" and choose this one. You can also click-and-drag your file in you tabular file reader.
