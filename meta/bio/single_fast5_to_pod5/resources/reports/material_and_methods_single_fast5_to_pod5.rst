Single fast5 to pod5 conversion
===============================

Samples which data was stored in a single archive of multiple single fast5 files were
downloaded with aria2c [#aria2c]_, extracted using GNU tar [#tar]_, aggregated to multi
fast5 files with ont_fast5_api [#ont_fast5_api]_, then converted to pod5 format with
pod5 [#pod5]_. Aria2c and tar are available on the snakemake [#snakemake]_ 
wrappers [#wrappers]_ through conda [#conda]_ environments, ont_fast5_api and pod5 
were deployed using apptainer [#apptainer]_ to solve undeclared dependencies.

.. [#aria2c] Tsujikawa Tatsuhiro et al. Aria2 [Computer software] https://github.com/aria2/aria2
.. [#tar] Gilmore John, et al. GNU tar Reference Manual: GNU tar an archiver tool. Samurai Media Limited, 2015.
.. [#ont_fast5_api] Brennen Forrest, et al. ont_fast5_api [computer software] https://github.com/nanoporetech/ont_fast5_api
.. [#pod5] Wang, Yunhao, et al. "Nanopore sequencing technology, bioinformatics and applications." Nature biotechnology 39.11 (2021): 1348-1365.
.. [#snakemake] Mölder, Felix, et al. "Sustainable data analysis with Snakemake." F1000Research 10 (2025): 33.
.. [#snakemake_wrappers] Köster, J, et al. Snakemake-Wrappers [Computer software] https://snakemake-wrappers.readthedocs.io/en/stable/
.. [#conda] Grüning, Björn, et al. "Bioconda: sustainable and comprehensive software distribution for the life sciences." Nature methods 15.7 (2018): 475-476.
.. [#apptainer] Pasupuleti, Rajesh, et al. "Apptainer-Based Containers for Legacy and Platform Dependent Machine Learning Applications on HPC Systems." 2024 IEEE International Conference on Big Data (BigData). IEEE, 2024. 
