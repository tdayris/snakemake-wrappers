Download genome sequence
========================

Ensembl [#ensembl]_ raw genome sequence was downloaded using Snakemake [#snakemake]_ 
wrapper's [#snakemake_wrappers]_ curl [#core_uutils]_, for {{ wildcards.species }},
build {{ wildcards.build }}, release {{ wildcards.release }}. Pyfaidx [#pyfaidx]_ 
was used to extract cannonical chromosomes. Samtools [#sambools]_, picard [#picard]_, 
and UCSC [#ucsc]_ were used to index and format genome sequence in commonly used 
formats accross bioinformatics. Xan [#xan]_ was used to extract chromosome sizes from 
fasta index.

.. [#ensembl] Yates, Andrew D., et al. "Ensembl 2026." Nucleic Acids Research 54.D1 (2026): D1053-D1060.
.. [#snakemake] Mölder, Felix, et al. "Sustainable data analysis with Snakemake." F1000Research 10 (2025): 33.
.. [#snakemake_wrappers] https://snakemake-wrappers.readthedocs.io/en/stable/
.. [#core_uutils] Ledru, Sylvestre, Samuel Tardieu, and Stefano Zacchiroli. "Rust Coreutils: Rebuilding Unix Foundations in a Modern Language." arXiv preprint arXiv:2608.07135 (2026).
.. [#pyfaidx] Shirley, Matthew D., et al. Efficient" pythonic" access to FASTA files using pyfaidx. No. e1196. PeerJ PrePrints, 2015.
.. [#sambools] Danecek, Petr, et al. "Twelve years of SAMtools and BCFtools." Gigascience 10.2 (2021): giab008.
.. [#picard] Van der Auwera, Geraldine A., and Brian D. O'Connor. Genomics in the cloud: using Docker, GATK, and WDL in Terra. O'Reilly Media, 2020.
.. [#ucsc] Navarro Gonzalez, Jairo, et al. "The UCSC genome browser database: 2021 update." Nucleic acids research 49.D1 (2021): D1046-D1057.
.. [#xan] Plique, Guillaume, et al. "xan, the CSV magician." (2026).
