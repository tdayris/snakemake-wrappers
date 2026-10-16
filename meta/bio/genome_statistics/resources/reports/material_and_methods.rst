Genome quality control
======================

Genome sequence and annotations qualities, formats and content were checked using
Snakemake [#snakemake]_ wrappers [#snakemake_wrappers]_. Agat [#agat]_, 
assembly-stats [#assembly_stats]_, seqkit [#seqkit]_, and Fastaguard [#fastaguard]_
were used sequencially on genome assemnbly. Results were formatted, and aggregated 
using go-yq [#go_yq]_ for easier parsing and accession. Quality report was built with 
MultiQC [#multiqc]_.

.. [snakemake] Mölder, Felix, et al. "Sustainable data analysis with Snakemake." F1000Research 10 (2025): 33.
.. [#snakemake_wrappers] Köster, J, et al. Snakemake-Wrappers [Computer software] https://snakemake-wrappers.readthedocs.io/en/stable/
.. [#agat] Dainat J. AGAT: Another Gff Analysis Toolkit to handle annotations in any GTF/GFF format. Zenodo. https://www.doi.org/10.5281/zenodo.3552717
.. [#assembly_stats] Thornton Jr, James. "Assembly Stats with MetaQuast." (2017).
.. [#seqkit] Shen, Wei, et al. "SeqKit: a cross-platform and ultrafast toolkit for FASTA/Q file manipulation." PloS one 11.10 (2016): e0163962.
.. [#fastaguard] ESTAJI, E. FastaGuard [Computer software]. https://github.com/ehsanestaji/FastaGuard
.. [#go_yq] Farah, M. YQ [Computer software]. https://github.com/mikefarah/yq

