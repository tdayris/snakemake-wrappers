Mapping with Dorado
===================

Dorado [#dorado]_ was used to call bases from pod5 [#pod5]_ and align reads from ubams.
Alignment was performed ensembl [#ensembl]_ on {{ wildcards.species }}, 
{{ wildcards.build }}, release {{ wildcards.release }}. Sambamba [#sambamba]_ was used 
to sort and index aligned reads. While sambamba was deployed through 
Snakemake [#snakemake]_ wrappers [#snakemake_wrappers]_ and pinned conda [#conda]_ 
environment, dorado and pod5 were deployed through apptainer [#apptainer]_ containers.

.. [#dorado] Oxford Nanopore PLC. Dorado [Computer software] https://github.com/nanoporetech/dorado
.. [#pod5] Wang, Yunhao, et al. "Nanopore sequencing technology, bioinformatics and applications." Nature biotechnology 39.11 (2021): 1348-1365.
.. [#ensembl] Yates, Andrew D., et al. "Ensembl 2026." Nucleic Acids Research 54.D1 (2026): D1053-D1060.
.. [#sambamba] Tarasov, Artem, et al. "Sambamba: fast processing of NGS alignment formats." Bioinformatics 31.12 (2015): 2032-2034.
.. [#snakemake] Mölder, Felix, et al. "Sustainable data analysis with Snakemake." F1000Research 10 (2025): 33.
.. [#snakemake_wrappers] Köster, J, et al. Snakemake-Wrappers [Computer software] https://snakemake-wrappers.readthedocs.io/en/stable/
.. [#conda] Grüning, Björn, et al. "Bioconda: sustainable and comprehensive software distribution for the life sciences." Nature methods 15.7 (2018): 475-476.
.. [#apptainer] Pasupuleti, Rajesh, et al. "Apptainer-Based Containers for Legacy and Platform Dependent Machine Learning Applications on HPC Systems." 2024 IEEE International Conference on Big Data (BigData). IEEE, 2024. 
