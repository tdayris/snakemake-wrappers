Fastq file accession
====================

Raw fastq files were accessed through Snakemake [#snakemake]_ storage plugins for
locally mounted filesystem [#plugin_fs]_ (using rsync [#rsync]_), distant files
accessible through HTTP [#plugin_http]_, and iRODS [#irods]_ protocol [#plugin_irods]_.

Samples that were splitted through multiple fastq were concatenated and checked using
fastcat [#fastcat]_, and compressed with crabz [#crabz]_. Sanity controls were made
using fastcat. Quality tables were aggregated and summarized using XSV [#xsv]_.

.. [#snakemake] Mölder, Felix, et al. "Sustainable data analysis with Snakemake." F1000Research 10 (2025): 33.
.. [#snakemake_wrappers] https://snakemake-wrappers.readthedocs.io/en/stable/
.. [#plugin_fs] Koester, J. snakemake-storage-plugin-fs [Computer software]. https://github.com/snakemake/snakemake-storage-plugin-fs
.. [#rsync] Trigdell, A. rsync [Computer software]. https://github.com/RsyncProject/rsync/
.. [#plugin_http] Tomkins-Tinch, C. snakemake-storage-plugin-http [Computer software]. https://github.com/snakemake/snakemake-storage-plugin-http
.. [#irods] Rajasekar, Arcot, et al. iRODS primer: integrated rule-oriented data system. Morgan & Claypool Publishers, 2010.
.. [#plugin_irods] Stolpe, O. snakemake-storage-plugin-irods [Computer software]. https://github.com/snakemake/snakemake-storage-plugin-irods
.. [#fastcat] Epi2me nanopore tech, fastcat [Computer software]. https://github.com/epi2me-labs/fastcat
.. [#crabz] Stadick, S. Crabz [Computer software]. https://github.com/sstadick/crabz
.. [#xsv] Gallant, A. XSV [Computer software]. https://github.com/BurntSushi/xsv
