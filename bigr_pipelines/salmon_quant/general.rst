Material and Methods:
#####################

Quality control were made on raw `FastP <https://en.wikipedia.org/wiki/FASTQ_format>`_ files with FastQC. Quality reports were gathered with MultiQC. Abundance estimation was performed with Salmon, which optional parameters (if any) are listed below.

Count aggregation and MultiQC customization were performed with in-house python scripts, which copies are available at the end of this report.

* Salmon index optional arguments: `{{snakemake.config.params.salmon_index_extra}}`
* Salmon quantification optional arguments: `{{snakemake.config.params.salmon_quant_extra}}`
* Library type was set to `{{snakemake.config.params.libType}}` in Salmon. See `Salmon's library type <https://salmon.readthedocs.io/en/latest/library_type.html>`_ definitions for more information.

The whole pipeline was powered by both `Snakemake <https://snakemake.readthedocs.io>`_ , and `Snakemake Wrappers <https://snakemake-wrappers.readthedocs.io/>`_ .

If you need any other information, please read the `Frequently Asked questions <https://github.com/tdayris-perso/rna-count-salmon#frequently-asked-questions-by-my-fellow-biologists-on-this-pipeline>`_ , then contact your bioinformatician if you're still in trouble.


Citations:
##########

This pipeline stands on best practices found in multiple high impact papers, published in Nature, Cell, Bioinformatics, and others.


FastQC
  ANDREWS, Simon, et al. FastQC: a quality control tool for high throughput sequence data. 2010.

  Why FastQC? FastQC is a quite popular tool among the field of bioinformatics and genomics. Cited more than 4000 times since its publication, this tool performs reliable quality controls on sequenced reads.

  https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

Salmon
  Patro, Rob, et al. “Salmon provides fast and bias-aware quantification of transcript expression.” Nature Methods (2017). Advanced Online Publication. doi: 10.1038/nmeth.4197.

  Why Salmon? Cited more than 1200 the last three years, Salmon is a very powerful tool to quantify transcripts taking into account a wide range of bias (either biological, or bioinformatical ones).

  https://salmon.readthedocs.io/

MultiQC
  EWELS, Philip, MAGNUSSON, Måns, LUNDIN, Sverker, et al. MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics, 2016, vol. 32, no 19, p. 3047-3048.

  Why MultiQC? MultiQC is a very efficient tool when it comes to quality gathering. It has been cited more than 500 times in a very wide range of journals including Nature, Bioinformatics, Cell, etc.

  https://multiqc.info/

Snakemake
  Köster, Johannes and Rahmann, Sven. “Snakemake - A scalable bioinformatics workflow engine”. Bioinformatics 2012.

  Why Snakemake? Snakemake is a very popular workflow manager in data science and bioinformatics. It has about three new citations per week within the scopes of biology, medicine and bioinformatics.

  https://snakemake.readthedocs.io/
  https://snakemake-wrappers.readthedocs.io/
