#!/usr/bin/env python3
# coding: utf-8

"""Write down material and methods for this pipeline"""

from typing import Optional, List


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


def immunedeconv_tool_list(org: str, nb) -> List[str]:
    dbs = []
    if org.lower().startswith("grch38"):
        dbs = [
            "xCell",
            "quanTIseq",
            "EPIC",
            "MCP-counter",
            "CIBERSORT_ABS",
            "CIBERSORT",
        ]
    else:
        dbs = [
            "mMCP-counter",
            "CIBERSORT_ABS",
            "CIBERSORT",
        ]
    result = [
        f"{db}[{nb + idx}]" for idx, db in enumerate(dbs)
    ]
    qts = [
        f"{quote_dict[db]}[{nb + idx}]" for idx, db in enumerate(dbs)
    ]

    yield from zip(result, qts)

quote_dict = {
    "fastp": [
        "0.23.2", 
        "Shifu Chen, Yanqing Zhou, Yaru Chen, Jia Gu; fastp: an ultra-fast all-in-one FASTQ preprocessor, Bioinformatics, Volume 34, Issue 17, 1 September 2018, Pages i884–i890, https://doi.org/10.1093/bioinformatics/bty560"
    ],
    "fastqscreen": [
        "0.5.2", 
        "FastqScreen: Wingett SW, Andrews S. FastQ Screen: A tool for multi-genome mapping and quality control. F1000Res. 2018 Aug 24;7:1338. doi: 10.12688/f1000research.15931.2. PMID: 30254741; PMCID: PMC6124377."
    ],
    "multiqc": [
        "1.12", 
        "Philip Ewels, Måns Magnusson, Sverker Lundin and Max Käller MultiQC: Summarize analysis results for multiple tools and samples in a single report. Bioinformatics (2016) doi: 10.1093/bioinformatics/btw354 PMID: 27312411"
    ],
    "salmon": [
        "1.8.0",
        "Patro, Rob, et al. \"Salmon provides fast and bias-aware quantification of transcript expression.\" Nature methods 14.4 (2017): 417-419."
    ],
    "salmon_method": [
        "XXX",
        'Srivastava, Avi, et al. "Alignment and mapping methodology influence transcript abundance estimation." Genome biology 21.1 (2020): 1-29.'
    ],
    "tximport": [
        "1.14.0",
        'Love, Michael I., et al. "Tximeta: Reference sequence checksums for provenance identification in RNA-seq." PLoS computational biology 16.2 (2020): e1007664.'
    ],
    "swimming": [
        "XXX",
        'Love, Michael I., Charlotte Soneson, and Rob Patro. "Swimming downstream: statistical analysis of differential transcript usage following Salmon quantification." F1000Research 7 (2018).'
    ],
    "deseq": [
        "1.34.0",
        'Love, Michael, Simon Anders, and Wolfgang Huber. "Differential analysis of count data–the DESeq2 package." Genome Biol 15.550 (2014): 10-1186.'
    ],
    "pcaexplorer": [
        "2.18.0",
        'Marini, Federico, and Harald Binder. "pcaExplorer: an R/Bioconductor package for interacting with RNA-seq principal components." BMC bioinformatics 20.1 (2019): 1-8.'
    ],
    "bioinfo-kit": [
        "1.0.8",
        'Huang, Yu, et al. "Bioinfo-Kit: A Sharing Software Tool for Bioinformatics." Applied Mechanics and Materials. Vol. 472. Trans Tech Publications Ltd, 2014.'
    ],
    "enhanced_volcano": [
        "1.4.0",
        'Blighe, K., Sharmila Rana, and Myles Lewis. "EnhancedVolcano: Publication-ready volcano plots with enhanced colouring and labeling. R package version 1.2. 0." (2019): 9.'
    ],
    "clusterprofiler": [
        "4.2",
        'Wu, Tianzhi, et al. "clusterProfiler 4.0: A universal enrichment tool for interpreting omics data." The Innovation 2.3 (2021): 100141.'
    ],
    "MSigDB": [
        "XXX",
        'Liberzon, Arthur, et al. "Molecular signatures database (MSigDB) 3.0." Bioinformatics 27.12 (2011): 1739-1740.'
    ],
    "Reactome": [
        "XXX",
        'Fabregat, Antonio, et al. "The reactome pathway knowledgebase." Nucleic acids research 44.D1 (2016): D481-D487.'
    ],
    "WikiPathways": [
        "XXX",
        'Kutmon, Martina, et al. "WikiPathways: capturing the full diversity of pathway knowledge." Nucleic acids research 44.D1 (2016): D488-D494.'
    ],
    'microRNA Targets': [
        "XXX",
        'John, Bino, et al. "Human microRNA targets." PLoS biology 2.11 (2004): e363.'
    ],
    'miRDB': [
        "XXX", 
        'Chen, Yuhao, and Xiaowei Wang. "miRDB: an online database for prediction of functional microRNA targets." Nucleic acids research 48.D1 (2020): D127-D131.'
    ],
    "Gene Onthology": [
        "XXX", 
        'Gene Ontology Consortium. "The gene ontology resource: 20 years and still GOing strong." Nucleic acids research 47.D1 (2019): D330-D338.'
    ],
    "Human Phenotype Onthology": [
        "XXX", 
        'Köhler, Sebastian, et al. "The human phenotype ontology in 2017." Nucleic acids research 45.D1 (2017): D865-D876.'
    ],
    "ImmuneSigDB": [
        "XXX", 
        'Godec, Jernej, et al. "Compendium of immune signatures identifies conserved and species-specific biology in response to inflammation." Immunity 44.1 (2016): 194-206.'
    ],
    "Human Protein Atlas": [
        "XXX", 
        'Pontén, Fredrik, et al. "The Human Protein Atlas as a proteomic resource for biomarker discovery." Journal of internal medicine 270.5 (2011): 428-446.'
    ],
    "CORUM": [
        "XXX", 
        'Ruepp, Andreas, et al. "CORUM: the comprehensive resource of mammalian protein complexes." Nucleic acids research 36.suppl_1 (2007): D646-D650.'
    ],
    "SMART": [
        "XXX", 
        'Letunic, Ivica, Tobias Doerks, and Peer Bork. "SMART 6: recent updates and new developments." Nucleic acids research 37.suppl_1 (2009): D229-D232.'
    ],
    "Pfam": [
        "XXX", 
        'Finn, Robert D., et al. "Pfam: the protein families database." Nucleic acids research 42.D1 (2014): D222-D230.'
    ],
    "Interpro": [
        "XXX", 
        'Blum, Matthias, et al. "The InterPro protein families and domains database: 20 years on." Nucleic acids research 49.D1 (2021): D344-D354.'
    ],
    "STRING": [
        "XXX",
        'Mering, Christian von, et al. "STRING: a database of predicted functional associations between proteins." Nucleic acids research 31.1 (2003): 258-261.'
    ],
    "Disease Ontology": [
        "XXX",
        'Kibbe, Warren A., et al. "Disease Ontology 2015 update: an expanded and updated database of human diseases for linking biomedical knowledge through disease data." Nucleic acids research 43.D1 (2015): D1071-D1078.'
    ],
    "Network of Cancer Genes": [
        "XXX",
        'An, Omer, et al. "NCG 4.0: the network of cancer genes in the era of massive mutational screenings of cancer genomes." Database 2014 (2014).'
    ],
    "DisGeNET": [
        "XXX",
        'Piñero, Janet, et al. "The DisGeNET knowledge platform for disease genomics: 2019 update." Nucleic acids research 48.D1 (2020): D845-D855.'
    ],
    "snakemake": [
        "7.15.2",
        'Köster, Johannes, and Sven Rahmann. "Snakemake—a scalable bioinformatics workflow engine." Bioinformatics 28.19 (2012): 2520-2522.'
    ],
    "snakemake_pipeline": [
        "1.19.2",
        'Mölder, Felix, et al. "Sustainable data analysis with Snakemake." F1000Research 10 (2021).'
    ],
    "Immunedeconv": [
        "2.0.3",
        'Sturm, Gregor, Francesca Finotello, and Markus List. "Immunedeconv: an R package for unified access to computational methods for estimating immune cell fractions from bulk RNA-sequencing data." Bioinformatics for cancer immunotherapy: methods and protocols (2020): 223-232.'
    ],
    "xCell": [
        "XXX",
        'Aran, Dvir, Zicheng Hu, and Atul J. Butte. "xCell: digitally portraying the tissue cellular heterogeneity landscape." Genome biology 18 (2017): 1-14.'
    ],
    "quanTIseq": [
        "XXX",
        'Plattner, Christina, Francesca Finotello, and Dietmar Rieder. "Deconvoluting tumor-infiltrating immune cells from RNA-seq data using quanTIseq." Methods in enzymology. Vol. 636. Academic Press, 2020. 261-285.'
    ],
    "EPIC": [
        "XXX",
        'Racle, Julien, and David Gfeller. "EPIC: a tool to estimate the proportions of different cell types from bulk gene expression data." Bioinformatics for Cancer Immunotherapy: Methods and Protocols (2020): 233-248.'
    ],
    "MCP-counter": [
        "XXX",
        'Becht, Etienne, et al. "Estimating the population abundance of tissue-infiltrating immune and stromal cell populations using gene expression." Genome biology 17.1 (2016): 1-20.'
    ],
    "CIBERSORT": [
        "XXX",
        'Newman, Aaron M., et al. "Robust enumeration of cell subsets from tissue expression profiles." Nature methods 12.5 (2015): 453-457.'
    ],
    "CIBERSORT_ABS": [
        "XXX",
        'Chen, Binbin, et al. "Profiling tumor infiltrating immune cells with CIBERSORT." Cancer Systems Biology: Methods and Protocols (2018): 243-259.'
    ],
    "estimate": [
        "XXX",
        'Yoshihara, K., Shahmoradgoli, M., Martínez, E., Vegesna, R., Kim, H., Torres-Garcia, W., Treviño, V., Shen, H., Laird, P. W., Levine, D. A., Carter, S. L., Getz, G., Stemke-Hale, K., Mills, G. B., & Verhaak, R. G. (2013). Inferring tumour purity and stromal and immune cell admixture from expression data. Nature communications, 4, 2612. https://doi.org/10.1038/ncomms3612'
    ],
    "timer": [
        "XXX",
        'Li, B., Severson, E., Pignon, J.-C., Zhao, H., Li, T., Novak, J., … Liu, X. S. (2016). Comprehensive analyses of tumor immunity: implications for cancer immunotherapy. Genome Biology, 17(1), 174. https://doi.org/10.1186/s13059-016-1028-7'
    ],
    "abis": [
        "XXX",
        'Monaco, G., Lee, B., Xu, W., Mustafah, S., Hwang, Y. Y., ..., Larbi, A. (2019). RNA-Seq Signatures Normalized by mRNA Abundance Allow Absolute Deconvolution of Human Immune Cell Types. Cell reports, 26(6), 1627–1640.e7. https://doi.org/10.1016/j.celrep.2019.01.041'
    ],
    "mMCP-counter": [
        "XXX",
        'Petitprez, Florent, et al. "The murine Microenvironment Cell Population counter method to estimate abundance of tissue-infiltrating immune and stromal cell populations in murine samples using gene expression." Genome medicine 12 (2020): 1-15.'
    ],
    "seqimucc": [
        "XXX",
        'Chen Z, Quan L, Huang A, Zhao Q, Yuan Y, Yuan X, Shen Q, Shang J, Ben Y, Qin FX, Wu A. seq-ImmuCC: Cell-Centric View of Tissue Transcriptome Measuring Cellular Compositions of Immune Microenvironment From Mouse RNA-Seq Data. Front Immunol. 2018 Jun 5;9:1286. doi: 10.3389/fimmu.2018.01286. PMID: 29922297; PMCID: PMC5996037.'
    ],
    "dcq": [
        "XXX",
        'Altboum Z, Steuerman Y, David E, Barnett-Itzhaki Z, Valadarsky L, Keren-Shaul H, et al. Digital cell quantification identifies global immune cell dynamics during influenza infection. Mol Syst Biol. 2014;10(2):720.'
    ],
    "base": [
        "XXX",
        'Varn, F. S., Andrews, E. H., Mullins, D. W., & Cheng, C. (2016). Integrative analysis of breast cancer reveals prognostic haematopoietic activity and patient-specific immune response profiles. Nature communications, 7, 10248. https://doi.org/10.1038/ncomms10248'
    ],
}

quote_number = 0
config = snakemake.params["all_parameters"]
text = [
    "# Material and Methods",
    "## Quality controls",
]
quotes = [
    "# Citations"
]

quote_number += 1
text.append(
    f"Fastq file trimming and raw quality control was performed by "
    f"Fastp[{quote_number}] (version: {quote_dict['fastp'][0]}), "
    f"using {fastp_parameters(config['fastp']['adapters'])}."
)
quotes.append(
    f"[{quote_number}]: {quote_dict['fastp'][1]}"
)

if config["step"].get("qc", False):
    quote_number += 1
    text.append(
        f"FastqScreen[{quote_number}] (version: {quote_dict['fastqscreen'][0]}) was "
        "ran through trimmed fastq files to assess their specificity to exptected genomes "
        f'(here: {config["reference"]["genome_name"]}).'
    )
    quotes.append(
        f"[{quote_number}]: {quote_dict['fastqscreen'][1]}"
    )

quote_number += 1
text.append(
    "Through the whole pipeline, quality control reports were aggregated in a single HTML "
    f"file using MultiQC[{quote_number}] (version: {quote_dict['multiqc'][0]}) using default "
    "parameters, and an in-house header. Tools and graphs that are not available natively "
    "on MultiQC were included through configuration file and custom graphs method, "
    "using in-house scripts."
)
quotes.append(
    f"[{quote_number}]: {quote_dict['multiqc'][1]}"
)

quote_number += 1
text.append("## Gene Expression")
text.append(
    f"Transcript abundance estimation was done with Salmon[{quote_number}] "
    f"(version: {quote_dict['salmon'][0]}), "
    f"using {parameters_explained(config['salmon']['salmon_quant_extra'])}. "
)
quotes.append(
    f"[{quote_number}]: {quote_dict['salmon'][1]}"
)

quote_number += 1
text.append(
    f"This estimation was done over Ensembl {config['ref']['genome']} "
    f"using decoy sequences[{quote_number}]. "
)
quotes.append(
    f"[{quote_number}]: {quote_dict['salmon_method'][1]}"
)

quote_number += 1
text.append(
    f"Transcript to gene expression has been performed by tximport[{quote_number}] "
    f"(version: {quote_dict['tximport']}), "
    f"using {parameters_explained(config['deseq2']['tximport_extra'])}, "
)
quotes.append(
    f"[{quote_number}]: {quote_dict['tximport'][1]}"
)

quote_number += 1
text.append(
    f"This whole method was described in Swimming Downstream protocol[{quote_number}]."
)
quotes.append(
    f"[{quote_number}]: {quote_dict['swimming'][1]}"
)

if config["step"].get("dge") or config["step"].get("gsea"):

    quote_number += 1
    text.append(
        f"Differential gene expression was performed with DESeq2[{quote_number}] "
        f"(version: {quote_dict['deseq'][0]}), using an "
        f"alpha threshold of {config['deseq2']['thresholds']['alpha']}. "
    )
    quotes.append(
        f"[{quote_number}]: {quote_dict['deseq'][1]}"
    )

    quote_number += 1
    tmp = f"DESeq2 results were described using pcaExplorer[{quote_number}] "
    tmp += f"(version {quote_dict['pcaexplorer'][0]}), "
    quotes.append(
        f"[{quote_number}]: {quote_dict['pcaexplorer'][1]}"
    )

    quote_number += 1
    tmp += f'bioinfokit[{quote_number}] (version: {quote_dict["bioinfo-kit"][0]}), '
    quotes.append(
        f'[{quote_number}]: {quote_dict["bioinfo-kit"][1]}'
    )

    quote_number += 1
    tmp += f'and EnhancedVolcano[{quote_number}] (version: {quote_dict["enhanced_volcano"][0]}) '
    quotes.append(
        f'[{quote_number}]: {quote_dict["enhanced_volcano"][1]}'
    )
    tmp += "through in-house scripts."

if config["step"].get("gsea"):
    text.append([tmp, "## Gene Set Enrichment Analysis"])
    quote_number += 1
    text.append(
        f'Gene set enrichment analysis was performed with ClusterProfiler[{quote_number}] '
        f'(version: {quote_dict["clusterprofiler"][0]}) using in-house scripts '
    )
    quotes.append(
        f'[{quote_number}]: {quote_dict["clusterprofiler"][1]}'
    )

    tmp = "and the following databases: "
    dbs = [
        "MSigDB",
        "Reactome",
        "WikiPathways",
        "microRNA Targets",
        "miRDB",
        "Gene Onthology",
        "Human Phenotype Onthology", 
        "ImmuneSigDB", 
        "Human Protein Atlas", 
        "CORUM", 
        "SMART", 
        "Pfam", 
        "Interpro", 
        "STRING", 
        "Disease Ontology", 
        "Network of Cancer Genes", 
        "DisGeNET"
    ]

    for idx, db in enumerate(dbs):
        quote_number += 1
        if idx == len(dbs) - 1:
            tmp = f"and {db}[{quote_number}]."
        else:
            tmp += f"{db}[{quote_number}], "
        quotes.append(
            f'[{quote_number}]: {quote_dict[db][1]}'
        )


if config["step"].get("immunedeconv"):

    quote_number += 1
    tmp = f'Cell type deconvolution has been done with ImmuneDeconv[{quote_number}]'
    tmp += f"(version: {quote_dict['Immunedeconv'][0]}) and the following list of methods: "
    quotes.append(f"[{quote_number}]: {quote_dict['Immunedeconv'][1]}")

    quote_number += 1
    for db, quote in immunedeconv_tool_list(config["organism"], quote_number):
        tmp += f'{quote_dict[db]}, '
        quotes.append(quote)
        quote_number += 1

text += [tmp, "## Pipeline"]

quote_number += 1
tmp = f'Ths whole pipeline was prowered with Snakemake[{quote_number}] '
tmp += f'(version: {quote_dict["snakemake"][0]}), and'
quotes.append(f"[{quote_number}]: {quote_dict['snakemake'][1]}")

quote_number += 1
tmp += f'Snakemake Wrappers[{quote_number}]'
tmp += f' (version: {quote_dict["snakemake_pipeline"][0]}).'
quotes.append(f"[{quote_number}]: {quote_dict['snakemake_pipeline'][1]}")

text.append(tmp)
text += quotes

with open(snakemake.output[0], "w") as mnm_stream:
    mnm_stream.write("\n".join(text))
