This is the result of the abundance estimation of sample {{ snakemake.wildcards.sample }}.

It contains five columns:

+-------------------+--------------------------------------------------+
| Column            | Usage                                            |
+===================+==================================================+
| Name              | The transcript ID                                |
+-------------------+--------------------------------------------------+
| Length            | The transcript length in the sequence annotation |
+-------------------+--------------------------------------------------+
| Effective Length  | The actual length within the RNASeq              |
+-------------------+--------------------------------------------------+
| NumReads          | The raw reads counts (do not use that value)     |
+-------------------+--------------------------------------------------+
| TPM               | The TPM normalized read counts                   |
+-------------------+--------------------------------------------------+

This is a TSV-formatted file. It can be opened in your favorite tabular file reader, like LibreOffice Calc, Excel, etc. Open your tabular file reader, then hit "open file" and choose this one. You can also click-and-drag your file in you tabular file reader.
