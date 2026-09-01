Download genome annotation
==========================

Ensembl [ensembl]_ raw genome sequence is was downloaded using Snakemake [snakemake]_ 
wrapper's [snakemake_wrappers]_ curl [core_uutils]_ for {{ wildcards.species }}, 
build {{ wildcards.build }}, release {{ wildcards.release }}. Agat [agat]_ 
configuration was build and tweaked using both Agat and go-yq [go_yq]_. Agat was used 
to remove genic object with a TSL below 5, fix common ensembl GTF/GFF format issue, and
ensure all contigs present in the genome sequence match this genome annotation.

.. [ensembl] Yates, Andrew D., et al. "Ensembl 2026." Nucleic Acids Research 54.D1 (2026): D1053-D1060.
.. [snakemake] Mölder, Felix, et al. "Sustainable data analysis with Snakemake." F1000Research 10 (2025): 33.
.. [snakemake_wrappers] https://snakemake-wrappers.readthedocs.io/en/stable/
.. [core_uutils] Ledru, Sylvestre, Samuel Tardieu, and Stefano Zacchiroli. "Rust Coreutils: Rebuilding Unix Foundations in a Modern Language." arXiv preprint arXiv:2608.07135 (2026).
.. [agat] Dainat J. AGAT: Another Gff Analysis Toolkit to handle annotations in any GTF/GFF format. Zenodo. https://www.doi.org/10.5281/zenodo.3552717
.. [go_yq] https://mikefarah.gitbook.io/yq
