#!/usr/bin/env python3
# coding: utf-8

"""Write down material and methods for this pipeline"""

from typing import Optional


def is_default(params: Optional[str] = None) -> bool:
    if params is None:
        return True
    elif params.lower() in ["none", "null", "None", "NULL", ""]:
        return True

    return False


def parameters_explained(params: Optional[str] = None) -> str:
    if is_default(params) is True:
        return "default arguments"
    return "the following default parameters: `{params}`"


def fastp_parameters(
    adapters: Optional[str] = None, extra: Optional[str] = None
) -> str:
    result = ""
    if not is_default(adapters):
        result = "custom adapter list, and using "
    result += parameters_explained(extra)

    return result


config = snakemake.params["all_parameters"]
text = f"""
# Material and Methods

## Quality controls

Fastq file trimming and raw quality control was performed by Fastp[1] (version: 0.20), using {fastp_parameters(config["fastp"]["adapters"])}. 
FastqScreen[2] (version: 0.5.2) was ran through trimmed fastq files to assess their specificity to exptected genomes (here: {config["reference"]["genome_name"]}).
Through the whole pipeline, quality control reports were aggregated in a single HTML file using MultiQC[3] (version: 1.10.1) using default parameters, and an in-house header. Tools and graphs that are not available natively on MultiQC were included through configuration file and custom graphs method, using in-house scripts.

## DGE

Transcript abundance estimation was done with Salmon[4] (version: 1.4.0), using {parameters_explained(config["salmon"]["salmon_quant_extra"])}. This estimation was done over Ensembl GRCh38.99 using decoy sequences[5]. Transcript to gene expression has been performed by tximport[6] (version: ), using {parameters_explained(config["deseq2"]["tximport_extra"])}, as described in [7].
Differential gene expression was performed with DESeq2[8] (version: 1.34.0), using an alpha threshold of {config["deseq2"]["thresholds"]["alpha"]}. DESeq2 results were described using pcaExplorer[9] (version 2.18.0), bioinfokit[10] (version: 1.0.8), and EnhancedVolcano[11] (version: 1.4.0) through in-house scripts.

## GSEA

Gene set enrichment analysis was performed with ClusterProfiler[12] (version: 4.2.0) using in-house scripts and the following databases: MSigDB[13], Reactome[14], WikiPathways[15], microRNA Targets[16], miRDB[17], Gene Onthology[18], Human Phenotype Onthology[19], ImmuneSigDB[20], Human Protein Atlas[21], CORUM[22], SMART[23], Pfam[24], Interpro[25], STRING[26], Disease Ontology[27], Network of Cancer Genes[28], and DisGeNET[29].

## Fusions

## MSI

## Deconvolution


## Pipeline

Ths whole pipeline was prowered with Snakemake[30] (version: 7.8.0), and Snakemake Wrappers[31] (version: 1.13.0).

# Citations

01. Shifu Chen, Yanqing Zhou, Yaru Chen, Jia Gu; fastp: an ultra-fast all-in-one FASTQ preprocessor, Bioinformatics, Volume 34, Issue 17, 1 September 2018, Pages i884–i890, https://doi.org/10.1093/bioinformatics/bty560
02. FastqScreen: Wingett SW, Andrews S. FastQ Screen: A tool for multi-genome mapping and quality control. F1000Res. 2018 Aug 24;7:1338. doi: 10.12688/f1000research.15931.2. PMID: 30254741; PMCID: PMC6124377.
03. Philip Ewels, Måns Magnusson, Sverker Lundin and Max Käller MultiQC: Summarize analysis results for multiple tools and samples in a single report. Bioinformatics (2016) doi: 10.1093/bioinformatics/btw354 PMID: 27312411 
04. Patro, Rob, et al. "Salmon provides fast and bias-aware quantification of transcript expression." Nature methods 14.4 (2017): 417-419.
05. Srivastava, Avi, et al. "Alignment and mapping methodology influence transcript abundance estimation." Genome biology 21.1 (2020): 1-29.
06. Love, Michael I., et al. "Tximeta: Reference sequence checksums for provenance identification in RNA-seq." PLoS computational biology 16.2 (2020): e1007664.
07. Love, Michael I., Charlotte Soneson, and Rob Patro. "Swimming downstream: statistical analysis of differential transcript usage following Salmon quantification." F1000Research 7 (2018).
08. Love, Michael, Simon Anders, and Wolfgang Huber. "Differential analysis of count data–the DESeq2 package." Genome Biol 15.550 (2014): 10-1186.
09. Marini, Federico, and Harald Binder. "pcaExplorer: an R/Bioconductor package for interacting with RNA-seq principal components." BMC bioinformatics 20.1 (2019): 1-8.
10. Huang, Yu, et al. "Bioinfo-Kit: A Sharing Software Tool for Bioinformatics." Applied Mechanics and Materials. Vol. 472. Trans Tech Publications Ltd, 2014.
11. Blighe, K., Sharmila Rana, and Myles Lewis. "EnhancedVolcano: Publication-ready volcano plots with enhanced colouring and labeling. R package version 1.2. 0." (2019): 9.
12. Wu, Tianzhi, et al. "clusterProfiler 4.0: A universal enrichment tool for interpreting omics data." The Innovation 2.3 (2021): 100141.
13. Liberzon, Arthur, et al. "Molecular signatures database (MSigDB) 3.0." Bioinformatics 27.12 (2011): 1739-1740.
14. Fabregat, Antonio, et al. "The reactome pathway knowledgebase." Nucleic acids research 44.D1 (2016): D481-D487.
15. Kutmon, Martina, et al. "WikiPathways: capturing the full diversity of pathway knowledge." Nucleic acids research 44.D1 (2016): D488-D494.
16. John, Bino, et al. "Human microRNA targets." PLoS biology 2.11 (2004): e363.
17. Chen, Yuhao, and Xiaowei Wang. "miRDB: an online database for prediction of functional microRNA targets." Nucleic acids research 48.D1 (2020): D127-D131.
18. Gene Ontology Consortium. "The gene ontology resource: 20 years and still GOing strong." Nucleic acids research 47.D1 (2019): D330-D338.
19. Köhler, Sebastian, et al. "The human phenotype ontology in 2017." Nucleic acids research 45.D1 (2017): D865-D876.
20. Godec, Jernej, et al. "Compendium of immune signatures identifies conserved and species-specific biology in response to inflammation." Immunity 44.1 (2016): 194-206.
21. Pontén, Fredrik, et al. "The Human Protein Atlas as a proteomic resource for biomarker discovery." Journal of internal medicine 270.5 (2011): 428-446.
22. Ruepp, Andreas, et al. "CORUM: the comprehensive resource of mammalian protein complexes." Nucleic acids research 36.suppl_1 (2007): D646-D650.
23. Letunic, Ivica, Tobias Doerks, and Peer Bork. "SMART 6: recent updates and new developments." Nucleic acids research 37.suppl_1 (2009): D229-D232.
24. Finn, Robert D., et al. "Pfam: the protein families database." Nucleic acids research 42.D1 (2014): D222-D230.
25. Blum, Matthias, et al. "The InterPro protein families and domains database: 20 years on." Nucleic acids research 49.D1 (2021): D344-D354.
26. Mering, Christian von, et al. "STRING: a database of predicted functional associations between proteins." Nucleic acids research 31.1 (2003): 258-261.
27. Kibbe, Warren A., et al. "Disease Ontology 2015 update: an expanded and updated database of human diseases for linking biomedical knowledge through disease data." Nucleic acids research 43.D1 (2015): D1071-D1078.
28. An, Omer, et al. "NCG 4.0: the network of cancer genes in the era of massive mutational screenings of cancer genomes." Database 2014 (2014).
29. Piñero, Janet, et al. "The DisGeNET knowledge platform for disease genomics: 2019 update." Nucleic acids research 48.D1 (2020): D845-D855.
30. Köster, Johannes, and Sven Rahmann. "Snakemake—a scalable bioinformatics workflow engine." Bioinformatics 28.19 (2012): 2520-2522.
31. Mölder, Felix, et al. "Sustainable data analysis with Snakemake." F1000Research 10 (2021).
"""

with open(snakemake.output[0], "w") as mnm_stream:
    mnm_stream.write(text)
