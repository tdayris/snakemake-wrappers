#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from collections import defaultdict
from json import dumps
from os import getcwd
from os.path import dirname, realpath

import logging

logging.basicConfig(
    filename=snakemake.log[0],
    filemode="w",
    level=logging.DEBUG
)

config = snakemake.params.get("config", defaultdict(str))
rundir = snakemake.params.get("rundir", f"data_{snakemake.wildcards['sample']}")
source = dirname(realpath(__file__))

sample = {
    "AnnotSV": {
        "DATABASE": "/mnt/beegfs/database/bioinfo/seqoia/original/annotations/AnnotSV/",
        "OPTIONS": "-genomeBuild GRCh38 -tx ENSEMBL",
        "VERSION": "2.5"
    },
    "CNV_FUSION_modificator": {
        "curie": "/mnt/beegfs/database/bioinfo/seqoia/original/scratch2/tmp/shared_files_tmp/CNV/gene_list/05-06-2018_WES_CancerGenesList_Large.txt",
        "gene_census": "/mnt/beegfs/database/bioinfo/seqoia/original/scratch2/tmp/shared_files_tmp/CNV/gene_list/cancer_gene_census_list_annotated.txt",
        "mercury": "/mnt/beegfs/database/bioinfo/seqoia/original/scratch2/tmp/shared_files_tmp/CNV/gene_list/cancer_gene_mercury.list.txt"
    },
    "Compte_rendu": {
        "QC": {
            #"MR363_WES-T": {
            f"{snakemake.wildcards['sample']}_WES-T": {
                "COUVERTURE": "-",
                "MEAN_COV": "223.5105312057919",
                "Q30": "-"
            },
            #"MR363_WGS-C": {
            f"{snakemake.wildcards['sample']}_WGS_C": {
                "COUVERTURE": "3.869986687038213",
                "MEAN_COV": "5.951654194084284",
                "Q30": "21708746826.0"
            },
            #"MR363_WTS": {
            f"{snakemake.wildcards['sample']}_WTS": {
                "COUVERTURE": "-",
                "MEAN_COV": "-",
                "Q30": "6465091868.0"
            },
            f"{snakemake.wildcards['sample']}_WGS-T": {
                "COUVERTURE": "-",
                "MEAN_COV": "-",
                "Q30": "-"
            }
        },
        "compte_rendu_content": "Versions: pipeline_cancer_wgs v2.1.0; snakefile_analysis_2.1; pipeline_config, v.2.1.0; cluster_config,v.2.1.0.\nMethode: Les fichiers bruts de sequencage (.BCL) du genome constitutionnel (WGS-C), de l exome tumoral (WES-T), du genome tumoral (WGS-T), et du transcriptome tumoral (WTS-T), si disponibles, sont generes grace a une etape de demultiplexage (bcl2fastq, v2.20.0.422, Illumina). Les sequences (format FASTQ) sont alignees contre le genome humain de reference (GRCh38.92.fa (GRCh38, release-92, Jul 02 2018, ftp.ensembl.org). Cette etape d alignement utilise une transformee de Burrows-Wheeler (BWA-MEM, 0.7.15). Un procede de nettoyage des fichiers d alignement est ensuite execute; celui-ci integre le marquage des duplicats de PCR (Picard MarkDuplicates (Picard Tools, 2.8.1)) et la recalibration des scores de qualite des bases (BaseRecalibrator, GATK4 (v4.1.0.0)).\nL appel des variants (SNP et Indels) sur le WGS-C et le WTS-T est realise par Haplotype Caller (GATK4, v4.1.0.0). L appel des variants (SNP et Indels) sur le WES-T est effectue par Mutect 2 (GATK4, v4.1.2.0).\nLes variants sont annotes au moyen de SNPeff(4.3t) et SnpSift(4.3t); les bases de donnees interrogees sont les suivantes : SNPEff (v4.3t), 1000Genomes (phase3, v2013-05-02), gnomAD exomes (v2.1.1), gnomAD genomes (v3), ClinVar (v20190722), COSMIC (coding, v89), COSMIC (non-coding, v89), dbscSNV (v1.1), dbSNP (v20180418), dbNSFP (v4.0), phastCons (v08-May-2015).\nLes variants somatiques sont utilises pour calculer la charge mutationnelle (pyTMB (v1.2.0)) et extraire les signatures mutationnelles (bases de donnees COSMIC (v3) integree a SigProfiler (v1.0.9)).\nLes fichiers d alignement du WGS-C et WGS-T sont exploites pour evaluer l instabilite des microsatellites (MSIsensor2 (v20191121)).\nLes CNVs sont detectes par Facet (v0.5.14) et WisecondorX (v1.1.5), puis annotes par AnnotSV (v2.5.1) auquel les bases de donnees suivantes ont ete ajoutees: Cytoband (Decembre 2013, USCS), COSMIC (v90).\nL appel des fusions est effectue independamment par Arriba (v1.2.0), STAR-Fusion (v1.9.0) et FusionCatcher (v1.10); les fusions sont validees par FusionInspector (STAR-Fusion v1.9.0) puis annotees par FusionAnnotator (STAR-Fusion, v1.9.0).",
        "pipeline_analyse_version": "SeqOIA_pipeline_cancer_wgs_v2.1.0",
        "pipeline_demux_version": "pipeline_demux_SeqOIA_v.3.2.2"
    },
    "QC": {
        "BED": "/mnt/beegfs/database/bioinfo/seqoia/original/annotations/Human/GRCh38/GRCh38.92.dict"
    },
    "QC_VALIDATION": {
        "WES": {
            "MEAN_COV": "150"
        },
        "WGS-C": {
            "COUVERTURE": "85",
            "MEAN_COV": "30",
            "Q30": "85000000000"
        },
        "WGS-T": {
            "COUVERTURE": "85",
            "MEAN_COV": "60",
            "Q30": "170000000000"
        },
        "WTS": {
            "Q30": "12800000000"
        }
    },
    "TMB": {
        "DB_YML": "/mnt/beegfs/database/bioinfo/seqoia/original/annotations/Human/GRCh38/DB/custom/TMB_curie/snpeff_jennifer_corrige.yml",
        "OPTIONS": "--effGenomeSize 37902905 --vaf 0.1 --maf 0.001 --minDepth 50 --minAltDepth 5 --filterLowQual --filterNonCoding --filterSyn",
        "VAR_YML": "/etc/curie-tmb/mutect2.yml"
    },
    "arriba_calling": {
        "DATABASE": "/database/blacklist_hg38_GRCh38_2018-11-04.tsv.gz",
        "GENOME_DIR": "/mnt/beegfs/database/bioinfo/seqoia/original/annotations/Human/GRCh38/index/star-fusion/v1.9/GRCh38_gencode_v33_CTAT_lib_Apr062020/ctat_genome_lib_build_dir/ref_annot.gtf.gz",
        "OPTIONS": "-T -P",
        "REFERENCE": "/mnt/beegfs/database/bioinfo/seqoia/original/annotations/Human/GRCh38/index/star-fusion/v1.9/GRCh38_gencode_v33_CTAT_lib_Apr062020/ctat_genome_lib_build_dir/ref_genome.fa.gz"
    },
    "arriba_drawing": {
        "CYTOBAND": "/database/cytobands_hg38_GRCh38_2018-02-23.tsv",
        "DATABASE": "/mnt/beegfs/database/bioinfo/seqoia/original/annotations/Human/GRCh38/index/star-fusion/v1.9/GRCh38_gencode_v33_CTAT_lib_Apr062020/ctat_genome_lib_build_dir/ref_annot.gtf.gz",
        "OPTIONS": "",
        "PROTEIN_ANNOT": "/database/protein_domains_hg38_GRCh38_2019-07-05.gff3"
    },
    "author": "Legendre Adrien",
    "bcl2fastq_wgs": {
        "OPTIONS": "--ignore-missing-positions --ignore-missing-controls --ignore-missing-filter --ignore-missing-bcls -r 4 -p 8 -w 4"
    },
    "bedtools": {
        "BED": "/mnt/beegfs/database/bioinfo/seqoia/original/annotations/Human/GRCh38/DB/custom/Twist_Bioscience/Extended_Boundaries/Merge_CORE-SPIKE_GRCh38.merged_bornes_10.bed",
        "OPTIONS": "-header -u"
    },
    "bwa_mem": {
        "OPTIONS": "-t 7 -M"
    },
    "description": "This config file is designed to run the SeqOIA pipeline cancer.",
    "facet_snp_pileup": {
        "DATABASE": "/mnt/beegfs/database/bioinfo/seqoia/original/annotations/Human/GRCh38/DB/public/dbsnp/20180418/All_20180418.vcf.gz",
        "OPTIONS": "-g -q15 -Q20 -P100 -r25,0"
    },
    "facet_snp_pileup_genome": {
        "DATABASE": "/mnt/beegfs/database/bioinfo/seqoia/original/annotations/Human/GRCh38/DB/public/dbsnp/20180418/All_20180418.vcf.gz",
        "OPTIONS": "-g -q15 -Q20 -P100 -r25,0"
    },
    "fastqc": {
        "OPTIONS": "--extract"
    },
    "fusion_annotator": {
        "GENOME_DIR": "/mnt/beegfs/database/bioinfo/seqoia/original/annotations/Human/GRCh38/index/star-fusion/v1.9/GRCh38_gencode_v33_CTAT_lib_Apr062020/ctat_genome_lib_build_dir/",
        "OPTIONS": ""
    },
    "fusion_catcher": {
        "ANNOT": "/mnt/beegfs/database/bioinfo/seqoia/original/annotations/Human/GRCh38/index/fusioncatcher/human_v90",
        "OPTIONS": "-p 16"
    },
    "gatk_applybqsr": {
        "OPTIONS": ""
    },
    "gatk_base_recalibrator": {
        "OPTIONS": ""
    },
    "gatk_callableloci": {
        "OPTIONS": "-mbq 10 -mmq 20 --maxDepth 10000 -mlmq 7 --minDepth 10"
    },
    "gatk_callableloci_wes": {
        "OPTIONS": "-mbq 10 -mmq 20 --maxDepth 10000 -mlmq 7 --minDepth 10"
    },
    "gatk_callableloci_wts": {
        "OPTIONS": "--filter_reads_with_N_cigar -mbq 10 -mmq 20 --maxDepth 10000 -mlmq 7 --minDepth 10"
    },
    "gatk_combinevariant": {
        "OPTIONS": "-genotypeMergeOptions UNIQUIFY"
    },
    "gatk_depth_coverage": {
        "OPTIONS": "--omitDepthOutputAtEachBase --omitIntervalStatistics -ct 1 -ct 20 -ct 30 -ct 60 -ct 150 -ct 200 -mmq 20 -mbq 10"
    },
    "gatk_depth_coverage_wes": {
        "OPTIONS": "--omitDepthOutputAtEachBase --omitIntervalStatistics -ct 1 -ct 20 -ct 30 -ct 60 -ct 150 -ct 200 -mmq 20 -mbq 10"
    },
    "gatk_depth_coverage_wes_gatk4": {
        "OPTIONS": "--omit-depth-output-at-each-base --omit-interval-statistics --summary-coverage-threshold 1 --summary-coverage-threshold 20 --summary-coverage-threshold 30 --summary-coverage-threshold 60 --summary-coverage-threshold 150 --summary-coverage-threshold 200 --read-filter MappingQualityReadFilter --minimum-mapping-quality 20 --min-base-quality 10"
    },
    "gatk_depth_coverage_wts": {
        "OPTIONS": "--filter_reads_with_N_cigar --omitDepthOutputAtEachBase --omitIntervalStatistics -ct 1 -ct 20 -ct 30 -ct 60 -ct 150 -ct 200 -mmq 20 -mbq 10"
    },
    "gatk_filter_mutect2_wes": {
        "OPTIONS": ""
    },
    "gatk_filter_mutect2_wgs": {
        "OPTIONS": ""
    },
    "gatk_genomicsDBImport": {
        "OPTIONS": ""
    },
    "gatk_genotype_gvcf": {
        "OPTIONS": "-stand-call-conf 30.0"
    },
    "gatk_haplotype_caller": {
        "OPTIONS": "-genotyping-mode DISCOVERY -sample-ploidy 2 -mbq 10 -stand-call-conf 30.0 --pcr-indel-model NONE"
    },
    "gatk_index": {
        "OPTIONS": ""
    },
    "gatk_mergevcf": {
        "OPTIONS": ""
    },
    "gatk_mutect2_wes": {
        "OPTIONS": "--smith-waterman FASTEST_AVAILABLE"
    },
    "gatk_mutect2_wes_tumor_only": {
        "OPTIONS": "--smith-waterman FASTEST_AVAILABLE"
    },
    "gatk_mutect2_wgs": {
        "OPTIONS": "--smith-waterman FASTEST_AVAILABLE"
    },
    "general_informations": {
        "ANALYSIS": "MR363",
        "BIOINFO_ANALYSIS_DATE": "15022022",
        "CHR": {
            "AUTOSOME": [
                "1",
                "2",
                "3",
                "4",
                "5",
                "6",
                "7",
                "8",
                "9",
                "10",
                "11",
                "12",
                "13",
                "14",
                "15",
                "16",
                "17",
                "18",
                "19",
                "20",
                "21",
                "22"
            ],
            "FACET": [
                "1",
                "2",
                "3",
                "4",
                "5",
                "6",
                "7",
                "8",
                "9",
                "10",
                "11",
                "12",
                "13",
                "14",
                "15",
                "16",
                "17",
                "18",
                "19",
                "20",
                "21",
                "22",
                "X"
            ],
            "IDENTITO_INDEX": [
                "1",
                "2",
                "3",
                "4",
                "5",
                "6",
                "7",
                "8",
                "9",
                "10",
                "11",
                "12",
                "13",
                "14",
                "15",
                "16",
                "17",
                "18",
                "19",
                "20",
                "22",
                "X",
                "Y"
            ],
            "INDEX": [
                "1",
                "2",
                "3",
                "4",
                "5",
                "6",
                "7",
                "8",
                "9",
                "10",
                "11",
                "12",
                "13",
                "14",
                "15",
                "16",
                "17",
                "18",
                "19",
                "20",
                "21",
                "22",
                "X",
                "Y",
                "MT"
            ],
            "INDEX_WES": [
                "1",
                "2",
                "3",
                "4",
                "5",
                "6",
                "7",
                "8",
                "9",
                "10",
                "11",
                "12",
                "13",
                "14",
                "15",
                "16",
                "17",
                "18",
                "19",
                "20",
                "21",
                "22",
                "X",
                "Y"
            ]
        },
        "CLUSTER_CONFIG": "/mnt/beegfs/scratch/bioinfo_core/B21081_LUFR_02/data_output/seqOIA/v.2.1.0.GRCh38.cluster.config.json",
        "FASTA_FILE": "/mnt/beegfs/database/bioinfo/seqoia/original/annotations/Human/GRCh38/index/sorted_primary_assemblies/GRCh38.92.fa",
        "ID_PIPELINE": "SeqOIA_pipeline_cancer_wgs_v2.1.0",
        "KNOWN_INDELS": "/mnt/beegfs/database/bioinfo/seqoia/original/annotations/Human/GRCh38/DB/custom/1000g/phase1/Mills_and_1000G_gold_standard.indels.hg38.filtered_nochr_final.vcf",
        "KNOWN_SNPS": "/mnt/beegfs/database/bioinfo/seqoia/original/annotations/Human/GRCh38/DB/custom/1000g/phase1/1000G_phase1.snps.high_confidence.hg38.filtered_nochr_final.vcf.gz",
        "PID_CRC": "0620A",
        "PIPELINE_CONFIG": str(snakemake.output["config"]),
        "PIPELINE_DEMUX": "pipeline_demux_SeqOIA_v.3.2.2",
        "PLATEFORME UNIT": "GCS SeqOIA",
        "REFERENCE_NAME": "GRCh38.92",
        "SAMPLES": [
            f"{snakemake.wildcards['sample']}_WGS-C_S1",
            f"{snakemake.wildcards['sample']}_WTS_S1",
            f"{snakemake.wildcards['sample']}_WGS-T_S1",
            f"{snakemake.wildcards['sample']}_WES-T_S1"
        ],
        "TARGET_TYPE": "Whole Genome",
        "VCF_ID": str(snakemake.wildcards["sample"]),
        "VCF_VERSION": "2.0",
        "WES-T": {
            "FLOWCELL": {
                "LANE": [
                    "1"
                ],
                "SURFACE": [
                    "1"
                ],
                "SWATH": [
                    "1"
                ]
            },
            #"ID": "MR363_WES-T_S1",
            "ID": f"{snakemake.wildcards['sample']}_WES-T_S1",
            "RUN": "22NNNN_ANNNNN_NNNN_AXXXXXXXXX",
            "SEQUENCEUR": "A00000"
        },
        "WGS-C": {
            "FLOWCELL": {
                "LANE": [
                    "1"
                ],
                "SURFACE": [
                    "1"
                ],
                "SWATH": [
                    "1"
                ]
            },
            #"ID": "MR363_WGS-C_S1",
            "ID": f"{snakemake.wildcards['sample']}_WGS-C_S1",
            "RUN": "22NNNN_ANNNNN_NNNN_AXXXXXXXXX",
            "SEQUENCEUR": "A00000"
        },
        "WGS-T": {
            "FLOWCELL": {
                "LANE": [
                    "1"
                ],
                "SURFACE": [
                    "1"
                ],
                "SWATH": [
                    "1"
                ]
            },
            #"ID": "MR363_WGS-T_S1",
            "ID": f"{snakemake.wildcards['sample']}_WGS-T_S1",
            "RUN": "22NNNN_ANNNNN_NNNN_AXXXXXXXXX",
            "SEQUENCEUR": "A00000"
        },
        "WTS": {
            "FLOWCELL": {
                "LANE": [
                    "1"
                ],
                "SURFACE": [
                    "1"
                ],
                "SWATH": [
                    "1"
                ]
            },
            #"ID": "MR363_WTS_S1",
            "ID": f"{snakemake.wildcards['sample']}_WTS_S1",
            "RUN": "22NNNN_ANNNNN_NNNN_AXXXXXXXXX",
            "SEQUENCEUR": "A00000"
        }
    },
    "general_path": {
        #"BASEDIR": "/mnt/beegfs/scratch/bioinfo_core/B21081_LUFR_02/data_output/seqOIA/data_MR363",
        "BASEDIR": rundir,
        #"CONF_PATH": "/mnt/beegfs/scratch/bioinfo_core/B21081_LUFR_02/data_output/seqOIA/",
        "CONF_PATH": getcwd(),
        #"DIR_PATH": "/mnt/beegfs/scratch/bioinfo_core/B21081_LUFR_02/data_output/seqOIA/data_MR363",
        # "DIR_PATH": rundir,
        "DIR_PATH": getcwd(),
        #"EXPORT_PATH": "/mnt/beegfs/scratch/bioinfo_core/B21081_LUFR_02/data_output/seqOIA/data_MR363",
        "EXPORT_PATH": rundir,
        #"INPUT_PATH": "/mnt/beegfs/scratch/bioinfo_core/B21081_LUFR_02/data_output/seqOIA/data_MR363",
        "INPUT_PATH": rundir,
        "IRODS_OUTPUT_PATH": "",
        #"SCRIPT_DIR": "/mnt/beegfs/scratch/bioinfo_core/B21081_LUFR_02/data_output/seqOIA/script",
        "SCRIPT_DIR": f"{source}/seqOIA_script",
        #"SNAKEMAKE_RULES": "/mnt/beegfs/scratch/bioinfo_core/B21081_LUFR_02/data_output/seqOIA/rules"
        "SNAKEMAKE_RULES": f"{source}/seqOIA_rules"
    },
    "identito_vigilance": {
        "OPTIONS": "-tg 99 -ts 97.5"
    },
    "msisensor": {
        "MICROSATELLITES_LIST": "/mnt/beegfs/database/bioinfo/seqoia/original/annotations/Human/GRCh38/index/msisensor/GRCh38.92.microstat",
        "OPTIONS": ""
    },
    "picard_collect_insert_size_metrics": {
        "OPTIONS": "--VALIDATION_STRINGENCY=LENIENT"
    },
    "picard_collect_pcr_metrics": {
        "OPTIONS": ""
    },
    "picard_collect_quality_yield_metrics": {
        "OPTIONS": "--USE_ORIGINAL_QUALITIES False"
    },
    "picard_collect_summary_metrics": {
        "OPTIONS": ""
    },
    "picard_collect_variant_calling_metrics": {
        "DATABASE": "/mnt/beegfs/database/bioinfo/seqoia/original/annotations/Human/GRCh38/DB/public/dbsnp/20180418/All_20180418.vcf.gz",
        "GENOME": "GRCh38",
        "OPTIONS": ""
    },
    "picard_collect_wgs_metrics": {
        "GENOME": "GRCh38",
        "OPTIONS": "-MQ 20 -Q 20 -CAP 10000"
    },
    "picard_markdup": {
        "OPTIONS": "REMOVE_DUPLICATES=false ASSUME_SORTED=true DUPLICATE_SCORING_STRATEGY=TOTAL_MAPPED_REFERENCE_LENGTH \\\nREAD_NAME_REGEX='[a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*.' OPTICAL_DUPLICATE_PIXEL_DISTANCE='100' \\\nVALIDATION_STRINGENCY=LENIENT QUIET=true VERBOSITY=ERROR"
    },
    "sambamba_flagstat": {
        "OPTIONS": "-t 8"
    },
    "sambamba_index": {
        "OPTIONS": "-t 6"
    },
    "sambamba_index_bam": {
        "OPTIONS": "-t 6"
    },
    "sambamba_index_genome": {
        "OPTIONS": "-t 8"
    },
    "sambamba_index_wgs_markdup": {
        "OPTIONS": "-t 6"
    },
    "sambamba_index_wts_markdup": {
        "OPTIONS": "-t 6"
    },
    "sambamba_merge": {
        "OPTIONS": "-t 8"
    },
    "sambamba_merge_facet": {
        "OPTIONS": "-t 8"
    },
    "sambamba_slice": {
        "OPTIONS": ""
    },
    "samtools_idxstats": {
        "OPTIONS": ""
    },
    "samtools_mpileup_wgs": {
        "OPTIONS": "-B -C 0 -d 1000000"
    },
    "samtools_sort": {
        "OPTIONS": "-m 2G"
    },
    "sign_mut_ext": {
        "OPTIONS": "-cpu 8"
    },
    "snpeff_tmb": {
        "DATABASE": "/mnt/beegfs/database/bioinfo/seqoia/original/annotations/Human/GRCh38/index/snpEff",
        "GENOME": "GRCh38.92",
        "OPTIONS": "-i vcf -o vcf -noInteraction",
        "VERSION": "v4.3t"
    },
    "snpeff_wgs": {
        "DATABASE": "/mnt/beegfs/database/bioinfo/seqoia/original/annotations/Human/GRCh38/index/snpEff",
        "GENOME": "GRCh38.92",
        "OPTIONS": "-i vcf -o vcf -noInteraction",
        "VERSION": "v4.3t"
    },
    "snpsift_annotate_leaves": {
        "DATABASES": [
            {
                "info_fields": "EAS_AF,EUR_AF,AFR_AF,AMR_AF,SAS_AF",
                "name": "kg_",
                "path": "/mnt/beegfs/database/bioinfo/seqoia/original/annotations/Human/GRCh38/DB/custom/1000g/phase3/1000GENOMES-phase_3_chr_{wildcards.index}.vcf.gz",
                "version": "phase3_20130502"
            },
            {
                "info_fields": "AC0,InbreedingCoeff,PASS,RF,allele_type,variant_type,was_mixed,n_alt_alleles,AC_raw,AC,AC_female,AC_male,AC_popmax,AC_afr,AC_afr_female,AC_afr_male,AC_amr,AC_amr_female,AC_amr_male,AC_asj,AC_asj_female,AC_asj_male,AC_eas,AC_eas_female,AC_eas_male,AC_fin,AC_fin_female,AC_fin_male,AC_nfe,AC_nfe_female,AC_nfe_male,AC_oth,AC_oth_female,AC_oth_male,AC_sas,AC_sas_female,AC_sas_male,AF_raw,AF,AF_female,AF_male,AF_popmax,popmax,AF_afr,AF_afr_female,AF_afr_male,AF_amr,AF_amr_female,AF_amr_male,AF_asj,AF_asj_female,AF_asj_male,AF_eas,AF_eas_female,AF_eas_male,AF_fin,AF_fin_female,AF_fin_male,AF_nfe,AF_nfe_male,AF_nfe_female,AF_oth,AF_oth_female,AF_oth_male,AF_sas,AF_sas_female,AF_sas_male,AN_raw,AN,AN_female,AN_male,AN_popmax,AN_afr,AN_afr_female,AN_afr_male,AN_amr,AN_amr_female,AN_amr_male,AN_asj,AN_asj_female,AN_asj_male,AN_eas,AN_eas_female,AN_eas_male,AN_fin,AN_fin_female,AN_fin_male,AN_nfe,AN_nfe_female,AN_nfe_male,AN_oth,AN_oth_female,AN_oth_male,AN_sas,AN_sas_female,AN_sas_male,controls_AC_raw,controls_AC,controls_AC_female,controls_AC_male,controls_AC_popmax,controls_AC_afr,controls_AC_afr_female,controls_AC_afr_male,controls_AC_amr,controls_AC_amr_female,controls_AC_amr_male,controls_AC_asj,controls_AC_asj_female,controls_AC_asj_male,controls_AC_eas,controls_AC_eas_female,controls_AC_eas_male,controls_AC_fin,controls_AC_fin_female,controls_AC_fin_male,controls_AC_nfe,controls_AC_nfe_female,controls_AC_nfe_male,controls_AC_oth,controls_AC_oth_female,controls_AC_oth_male,controls_AC_sas,controls_AC_sas_female,controls_AC_sas_male,controls_AF_raw,controls_AF,controls_AF_female,controls_AF_male,controls_AF_popmax,controls_popmax,controls_AF_afr,controls_AF_afr_female,controls_AF_afr_male,controls_AF_amr,controls_AF_amr_female,controls_AF_amr_male,controls_AF_asj,controls_AF_asj_female,controls_AF_asj_male,controls_AF_eas,controls_AF_eas_female,controls_AF_eas_male,controls_AF_fin,controls_AF_fin_female,controls_AF_fin_male,controls_AF_nfe,controls_AF_nfe_female,controls_AF_nfe_male,controls_AF_oth,controls_AF_oth_female,controls_AF_oth_male,controls_AF_sas,controls_AF_sas_female,controls_AF_sas_male,controls_AN_raw,controls_AN,controls_AN_female,controls_AN_male,controls_AN_popmax,controls_AN_afr,controls_AN_afr_female,controls_AN_afr_male,controls_AN_amr,controls_AN_amr_female,controls_AN_amr_male,controls_AN_asj,controls_AN_asj_female,controls_AN_asj_male,controls_AN_eas,controls_AN_eas_female,controls_AN_eas_male,controls_AN_fin,controls_AN_fin_female,controls_AN_fin_male,controls_AN_nfe,controls_AN_nfe_female,controls_AN_nfe_male,controls_AN_oth,controls_AN_oth_female,controls_AN_oth_male,controls_AN_sas,controls_AN_sas_female,controls_AN_sas_male,controls_nhomalt_raw,controls_nhomalt,controls_nhomalt_female,controls_nhomalt_male,controls_nhomalt_popmax,controls_nhomalt_afr,controls_nhomalt_afr_female,controls_nhomalt_afr_male,controls_nhomalt_amr,controls_nhomalt_amr_female,controls_nhomalt_amr_male,controls_nhomalt_asj,controls_nhomalt_asj_female,controls_nhomalt_asj_male,controls_nhomalt_eas,controls_nhomalt_eas_female,controls_nhomalt_eas_male,controls_nhomalt_fin,controls_nhomalt_fin_female,controls_nhomalt_fin_male,controls_nhomalt_nfe,controls_nhomalt_nfe_female,controls_nhomalt_nfe_male,controls_nhomalt_oth,controls_nhomalt_oth_female,controls_nhomalt_oth_male,controls_nhomalt_sas,controls_nhomalt_sas_female,controls_nhomalt_sas_male,nhomalt_raw,nhomalt,nhomalt_female,nhomalt_male,nhomalt_popmax,nhomalt_afr,nhomalt_afr_female,nhomalt_afr_male,nhomalt_amr,nhomalt_amr_female,nhomalt_amr_male,nhomalt_asj,nhomalt_asj_female,nhomalt_asj_male,nhomalt_eas,nhomalt_eas_female,nhomalt_eas_male,nhomalt_fin,nhomalt_fin_female,nhomalt_fin_male,nhomalt_nfe,nhomalt_nfe_female,nhomalt_nfe_male,nhomalt_oth,nhomalt_oth_female,nhomalt_oth_male,nhomalt_sas,nhomalt_sas_female,nhomalt_sas_male,non_cancer_AC_raw,non_cancer_AC,non_cancer_AC_female,non_cancer_AC_male,non_cancer_AC_popmax,non_cancer_AC_afr,non_cancer_AC_afr_female,non_cancer_AC_afr_male,non_cancer_AC_amr,non_cancer_AC_amr_female,non_cancer_AC_amr_male,non_cancer_AC_asj,non_cancer_AC_asj_female,non_cancer_AC_asj_male,non_cancer_AC_eas,non_cancer_AC_eas_female,non_cancer_AC_eas_male,non_cancer_AC_fin,non_cancer_AC_fin_female,non_cancer_AC_fin_male,non_cancer_AC_nfe,non_cancer_AC_nfe_female,non_cancer_AC_nfe_male,non_cancer_AC_oth,non_cancer_AC_oth_female,non_cancer_AC_oth_male,non_cancer_AC_sas,non_cancer_AC_sas_female,non_cancer_AC_sas_male,non_cancer_AF_raw,non_cancer_AF,non_cancer_AF_female,non_cancer_AF_male,non_cancer_AF_popmax,non_cancer_popmax,non_cancer_AF_afr,non_cancer_AF_afr_female,non_cancer_AF_afr_male,non_cancer_AF_amr,non_cancer_AF_amr_female,non_cancer_AF_amr_male,non_cancer_AF_asj,non_cancer_AF_asj_female,non_cancer_AF_asj_male,non_cancer_AF_eas,non_cancer_AF_eas_female,non_cancer_AF_eas_male,non_cancer_AF_fin,non_cancer_AF_fin_female,non_cancer_AF_fin_male,non_cancer_AF_nfe,non_cancer_AF_nfe_female,non_cancer_AF_nfe_male,non_cancer_AF_oth,non_cancer_AF_oth_female,non_cancer_AF_oth_male,non_cancer_AF_sas,non_cancer_AF_sas_female,non_cancer_AF_sas_male,non_cancer_AN_raw,non_cancer_AN,non_cancer_AN_female,non_cancer_AN_male,non_cancer_AN_popmax,non_cancer_AN_afr,non_cancer_AN_afr_female,non_cancer_AN_afr_male,non_cancer_AN_amr,non_cancer_AN_amr_female,non_cancer_AN_amr_male,non_cancer_AN_asj,non_cancer_AN_asj_female,non_cancer_AN_asj_male,non_cancer_AN_eas,non_cancer_AN_eas_female,non_cancer_AN_eas_male,non_cancer_AN_fin,non_cancer_AN_fin_female,non_cancer_AN_fin_male,non_cancer_AN_nfe,non_cancer_AN_nfe_female,non_cancer_AN_nfe_male,non_cancer_AN_oth,non_cancer_AN_oth_female,non_cancer_AN_oth_male,non_cancer_AN_sas,non_cancer_AN_sas_female,non_cancer_AN_sas_male,non_cancer_nhomalt_raw,non_cancer_nhomalt,non_cancer_nhomalt_female,non_cancer_nhomalt_male,non_cancer_nhomalt_popmax,non_cancer_nhomalt_afr,non_cancer_nhomalt_afr_female,non_cancer_nhomalt_afr_male,non_cancer_nhomalt_amr,non_cancer_nhomalt_amr_female,non_cancer_nhomalt_amr_male,non_cancer_nhomalt_asj,non_cancer_nhomalt_asj_female,non_cancer_nhomalt_asj_male,non_cancer_nhomalt_eas,non_cancer_nhomalt_eas_female,non_cancer_nhomalt_eas_male,non_cancer_nhomalt_fin,non_cancer_nhomalt_fin_female,non_cancer_nhomalt_fin_male,non_cancer_nhomalt_nfe,non_cancer_nhomalt_nfe_female,non_cancer_nhomalt_nfe_male,non_cancer_nhomalt_oth,non_cancer_nhomalt_oth_female,non_cancer_nhomalt_oth_male,non_cancer_nhomalt_sas,non_cancer_nhomalt_sas_female,non_cancer_nhomalt_sas_male,non_neuro_AC_raw,non_neuro_AC,non_neuro_AC_female,non_neuro_AC_male,non_neuro_AC_popmax,non_neuro_AC_afr,non_neuro_AC_afr_female,non_neuro_AC_afr_male,non_neuro_AC_amr,non_neuro_AC_amr_female,non_neuro_AC_amr_male,non_neuro_AC_asj,non_neuro_AC_asj_female,non_neuro_AC_asj_male,non_neuro_AC_eas,non_neuro_AC_eas_female,non_neuro_AC_eas_male,non_neuro_AC_fin,non_neuro_AC_fin_female,non_neuro_AC_fin_male,non_neuro_AC_nfe,non_neuro_AC_nfe_female,non_neuro_AC_nfe_male,non_neuro_AC_oth,non_neuro_AC_oth_female,non_neuro_AC_oth_male,non_neuro_AC_sas,non_neuro_AC_sas_female,non_neuro_AC_sas_male,non_neuro_AF_raw,non_neuro_AF,non_neuro_AF_female,non_neuro_AF_male,non_neuro_AF_popmax,non_neuro_popmax,non_neuro_AF_afr,non_neuro_AF_afr_female,non_neuro_AF_afr_male,non_neuro_AF_amr,non_neuro_AF_amr_female,non_neuro_AF_amr_male,non_neuro_AF_asj,non_neuro_AF_asj_female,non_neuro_AF_asj_male,non_neuro_AF_eas,non_neuro_AF_eas_female,non_neuro_AF_eas_male,non_neuro_AF_fin,non_neuro_AF_fin_female,non_neuro_AF_fin_male,non_neuro_AF_nfe,non_neuro_AF_nfe_female,non_neuro_AF_nfe_male,non_neuro_AF_oth,non_neuro_AF_oth_female,non_neuro_AF_oth_male,non_neuro_AF_sas,non_neuro_AF_sas_female,non_neuro_AF_sas_male,non_neuro_AN_raw,non_neuro_AN,non_neuro_AN_female,non_neuro_AN_male,non_neuro_AN_popmax,non_neuro_AN_afr,non_neuro_AN_afr_female,non_neuro_AN_afr_male,non_neuro_AN_amr,non_neuro_AN_amr_female,non_neuro_AN_amr_male,non_neuro_AN_asj,non_neuro_AN_asj_female,non_neuro_AN_asj_male,non_neuro_AN_eas,non_neuro_AN_eas_female,non_neuro_AN_eas_male,non_neuro_AN_fin,non_neuro_AN_fin_female,non_neuro_AN_fin_male,non_neuro_AN_nfe,non_neuro_AN_nfe_female,non_neuro_AN_nfe_male,non_neuro_AN_oth,non_neuro_AN_oth_female,non_neuro_AN_oth_male,non_neuro_AN_sas,non_neuro_AN_sas_female,non_neuro_AN_sas_male,non_neuro_nhomalt,non_neuro_nhomalt_female,non_neuro_nhomalt_male,non_neuro_nhomalt_popmax,non_neuro_nhomalt_raw,non_neuro_nhomalt_afr,non_neuro_nhomalt_afr_female,non_neuro_nhomalt_afr_male,non_neuro_nhomalt_amr,non_neuro_nhomalt_amr_female,non_neuro_nhomalt_amr_male,non_neuro_nhomalt_asj,non_neuro_nhomalt_asj_female,non_neuro_nhomalt_asj_male,non_neuro_nhomalt_eas,non_neuro_nhomalt_eas_female,non_neuro_nhomalt_eas_male,non_neuro_nhomalt_fin,non_neuro_nhomalt_fin_female,non_neuro_nhomalt_fin_male,non_neuro_nhomalt_nfe,non_neuro_nhomalt_nfe_female,non_neuro_nhomalt_nfe_male,non_neuro_nhomalt_oth,non_neuro_nhomalt_oth_female,non_neuro_nhomalt_oth_male,non_neuro_nhomalt_sas,non_neuro_nhomalt_sas_female,non_neuro_nhomalt_sas_male,non_topmed_AC,non_topmed_AC_female,non_topmed_AC_male,non_topmed_AC_popmax,non_topmed_AC_raw,non_topmed_AC_afr,non_topmed_AC_afr_female,non_topmed_AC_afr_male,non_topmed_AC_amr,non_topmed_AC_amr_female,non_topmed_AC_amr_male,non_topmed_AC_asj,non_topmed_AC_asj_female,non_topmed_AC_asj_male,non_topmed_AC_eas,non_topmed_AC_eas_female,non_topmed_AC_eas_male,non_topmed_AC_fin,non_topmed_AC_fin_female,non_topmed_AC_fin_male,non_topmed_AC_nfe,non_topmed_AC_nfe_female,non_topmed_AC_nfe_male,non_topmed_AC_oth,non_topmed_AC_oth_female,non_topmed_AC_oth_male,non_topmed_AC_sas,non_topmed_AC_sas_female,non_topmed_AC_sas_male,non_topmed_AF_raw,non_topmed_AF,non_topmed_AF_female,non_topmed_AF_male,non_topmed_AF_popmax,non_topmed_AF_afr,non_topmed_AF_afr_female,non_topmed_AF_afr_male,non_topmed_AF_amr,non_topmed_AF_amr_female,non_topmed_AF_amr_male,non_topmed_AF_asj,non_topmed_AF_asj_female,non_topmed_AF_asj_male,non_topmed_AF_eas,non_topmed_AF_eas_female,non_topmed_AF_eas_male,non_topmed_AF_fin,non_topmed_AF_fin_female,non_topmed_AF_fin_male,non_topmed_AF_nfe,non_topmed_AF_nfe_female,non_topmed_AF_nfe_male,non_topmed_AF_oth,non_topmed_AF_oth_female,non_topmed_AF_oth_male,non_topmed_AF_sas,non_topmed_AF_sas_female,non_topmed_AF_sas_male,non_topmed_AN_raw,non_topmed_AN,non_topmed_AN_female,non_topmed_AN_male,non_topmed_AN_popmax,non_topmed_popmax,non_topmed_AN_afr,non_topmed_AN_afr_female,non_topmed_AN_afr_male,non_topmed_AN_amr,non_topmed_AN_amr_female,non_topmed_AN_amr_male,non_topmed_AN_asj,non_topmed_AN_asj_female,non_topmed_AN_asj_male,non_topmed_AN_eas,non_topmed_AN_eas_female,non_topmed_AN_eas_male,non_topmed_AN_fin,non_topmed_AN_fin_female,non_topmed_AN_fin_male,non_topmed_AN_nfe,non_topmed_AN_nfe_female,non_topmed_AN_nfe_male,non_topmed_AN_oth,non_topmed_AN_oth_female,non_topmed_AN_oth_male,non_topmed_AN_sas,non_topmed_AN_sas_female,non_topmed_AN_sas_male,non_topmed_nhomalt_raw,non_topmed_nhomalt,non_topmed_nhomalt_female,non_topmed_nhomalt_male,non_topmed_nhomalt_popmax,non_topmed_nhomalt_afr,non_topmed_nhomalt_afr_female,non_topmed_nhomalt_afr_male,non_topmed_nhomalt_amr,non_topmed_nhomalt_amr_female,non_topmed_nhomalt_amr_male,non_topmed_nhomalt_asj,non_topmed_nhomalt_asj_female,non_topmed_nhomalt_asj_male,non_topmed_nhomalt_eas,non_topmed_nhomalt_eas_female,non_topmed_nhomalt_eas_male,non_topmed_nhomalt_fin,non_topmed_nhomalt_fin_female,non_topmed_nhomalt_fin_male,non_topmed_nhomalt_nfe,non_topmed_nhomalt_nfe_female,non_topmed_nhomalt_nfe_male,non_topmed_nhomalt_oth,non_topmed_nhomalt_oth_female,non_topmed_nhomalt_oth_male,non_topmed_nhomalt_sas,non_topmed_nhomalt_sas_female,non_topmed_nhomalt_sas_male,age_hist_het_bin_freq,age_hist_hom_bin_freq,dp_hist_all_bin_freq,dp_hist_alt_bin_freq",
                "name": "gnomAD_exomes_",
                "path": "/mnt/beegfs/database/bioinfo/seqoia/original/annotations/Human/GRCh38/DB/public/gnomad/2.1.1/gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz",
                "version": "2.1.1"
            },
            {
                "info_fields": "AC0,AS_VQSR,InbreedingCoeff,PASS,variant_type,n_alt_alleles,AC_raw,AN_raw,AF_raw,nhomalt_raw,AC,AN,AF,nhomalt,AC_female,AN_female,AF_female,nhomalt_female,AC_male,AN_male,AF_male,nhomalt_male,AC_asj,AN_asj,AF_asj,nhomalt_asj,AC_asj_female,AN_asj_female,AF_asj_female,nhomalt_asj_female,AC_asj_male,AN_asj_male,AF_asj_male,nhomalt_asj_male,AC_fin,AN_fin,AF_fin,nhomalt_fin,AC_fin_female,AN_fin_female,AF_fin_female,nhomalt_fin_female,AC_fin_male,AN_fin_male,AF_fin_male,nhomalt_fin_male,AC_ami,AN_ami,AF_ami,nhomalt_ami,AC_ami_female,AN_ami_female,AF_ami_female,nhomalt_ami_female,AC_ami_male,AN_ami_male,AF_ami_male,nhomalt_ami_male,AC_oth,AN_oth,AF_oth,nhomalt_oth,AC_oth_female,AN_oth_female,AF_oth_female,nhomalt_oth_female,AC_oth_male,AN_oth_male,AF_oth_male,nhomalt_oth_male,AC_afr,AN_afr,AF_afr,nhomalt_afr,AC_afr_female,AN_afr_female,AF_afr_female,nhomalt_afr_female,AC_afr_male,AN_afr_male,AF_afr_male,nhomalt_afr_male,AC_eas,AN_eas,AF_eas,nhomalt_eas,AC_eas_female,AN_eas_female,AF_eas_female,nhomalt_eas_female,AC_eas_male,AN_eas_male,AF_eas_male,nhomalt_eas_male,AC_sas,AN_sas,AF_sas,nhomalt_sas,AC_sas_female,AN_sas_female,AF_sas_female,nhomalt_sas_female,AC_sas_male,AN_sas_male,AF_sas_male,nhomalt_sas_male,AC_nfe,AN_nfe,AF_nfe,nhomalt_nfe,AC_nfe_female,AN_nfe_female,AF_nfe_female,nhomalt_nfe_female,AC_nfe_male,AN_nfe_male,AF_nfe_male,nhomalt_nfe_male,AC_amr,AN_amr,AF_amr,nhomalt_amr,AC_amr_female,AN_amr_female,AF_amr_female,nhomalt_amr_female,AC_amr_male,AN_amr_male,AF_amr_male,nhomalt_amr_male",
                "name": "gnomAD_genomes_",
                "path": "/mnt/beegfs/database/bioinfo/seqoia/original/annotations/Human/GRCh38/DB/custom/gnomad/3.0/gnomad.genomes.r3.0.sites.chr{wildcards.index}.vcf.gz",
                "version": "3.0"
            },
            {
                "info_fields": "ALLELEID,CLNSIG,CLNREVSTAT",
                "name": "clinvar_",
                "path": "/mnt/beegfs/database/bioinfo/seqoia/original/annotations/Human/GRCh38/DB/custom/clinvar/20190722/clinvar_20190722_chr_{wildcards.index}.vcf.gz",
                "version": "20190722"
            },
            {
                "info_fields": "RS",
                "name": "dbsnp_",
                "path": "/mnt/beegfs/database/bioinfo/seqoia/original/annotations/Human/GRCh38/DB/custom/dbsnp/20180418/All_20180418_chr_{wildcards.index}.vcf.gz",
                "version": "20180418"
            },
            {
                "info_fields": "ID",
                "name": "cosmic_coding_",
                "path": "/mnt/beegfs/database/bioinfo/seqoia/original/annotations/Human/GRCh38/DB/custom/cosmic/v89/CosmicCodingMuts.custom.gunzip_chr_{wildcards.index}.vcf.gz",
                "version": "v89"
            },
            {
                "info_fields": "ID",
                "name": "cosmic_noncoding_",
                "path": "/mnt/beegfs/database/bioinfo/seqoia/original/annotations/Human/GRCh38/DB/custom/cosmic/v89/CosmicNonCodingVariants.custom.gunzip_chr_{wildcards.index}.vcf.gz",
                "version": "v89"
            },
            {
                "info_fields": "ADA_SCORE,RF_SCORE",
                "name": "dbscsnv_",
                "path": "/mnt/beegfs/database/bioinfo/seqoia/original/annotations/Human/GRCh38/DB/custom/dbscsnv/1.1/hg38_dbscsnv11_chr_{wildcards.index}.vcf.gz",
                "version": "v1.1"
            }
        ],
        "OPTIONS": "-noDownload -noId -tabix"
    },
    "snpsift_annotate_tmb": {
        "DATABASES": "/mnt/beegfs/database/bioinfo/seqoia/original/annotations/Human/GRCh38/DB/custom/gnomad/3.0/gnomad.genomes.r3.0.sites.vcf.bgz",
        "INFOS": "AF",
        "NAME": "gnomAD_genomes_",
        "OPTIONS": "-noDownload -noId -tabix",
        "VERSION": "3.0"
    },
    "snpsift_annotate_wes": {
        "OPTIONS": "-noDownload -noId -tabix -name Mutect2_ -info AF"
    },
    "snpsift_dbnsfp": {
        "DATABASE": "/mnt/beegfs/database/bioinfo/seqoia/original/annotations/Human/hg38/DB/public/dbnsfp/4.0a/dbNSFP4.0a.txt.gz",
        "OPTIONS": "-v -collapse -f SIFT_score,SIFT_pred,Polyphen2_HDIV_score,Polyphen2_HDIV_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_pred,LRT_score,LRT_pred,MutationTaster_score,MutationTaster_pred,FATHMM_score,FATHMM_pred,MetaSVM_score,MetaSVM_pred,MetaLR_score,MetaLR_pred,PROVEAN_score,PROVEAN_pred,CADD_raw,CADD_phred,GERP++_NR,GERP++_RS"
    },
    "snpsift_phastcons": {
        "DATABASE": "/mnt/beegfs/database/bioinfo/seqoia/original/annotations/Human/hg38/index/phastCons/",
        "OPTIONS": ""
    },
    "star": {
        "GENOME_DIR": "/mnt/beegfs/database/bioinfo/seqoia/original/annotations/Human/GRCh38/index/sorted_primary_assemblies/star/2.7.2d/",
        "OPTIONS": "--runThreadN 8 --outSAMtype BAM SortedByCoordinate --runMode alignReads --outSJfilterReads Unique --chimSegmentMin 10 --outSAMmapqUnique 60 --readFilesCommand zcat --twopassMode Basic --twopass1readsN -1"
    },
    "star_arriba": {
        "GENOME_DIR": "/mnt/beegfs/database/bioinfo/seqoia/original/annotations/Human/GRCh38/index/star/seqoia/2.7.5a/",
        "OPTIONS": "--runThreadN 8 --genomeLoad NoSharedMemory --outStd BAM_Unsorted --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outBAMcompression 0 --outFilterMultimapNmax 1 --outFilterMismatchNmax 3 --chimSegmentMin 10 --chimOutType SeparateSAMold --chimJunctionOverhangMin 10 --chimScoreMin 1 --chimScoreDropMax 30 --chimScoreJunctionNonGTAG 0 --chimScoreSeparation 1 --alignSJstitchMismatchNmax 5 -1 5 5 --chimSegmentReadGapMax 3"
    },
    "star_fusion": {
        "GENOME_DIR": "/mnt/beegfs/database/bioinfo/seqoia/original/annotations/Human/GRCh38/index/star-fusion/v1.9/GRCh38_gencode_v33_CTAT_lib_Apr062020/ctat_genome_lib_build_dir/",
        "OPTIONS": "--CPU 12"
    },
    "star_fusion_inspector": {
        "GENOME_DIR": "/mnt/beegfs/database/bioinfo/seqoia/original/annotations/Human/GRCh38/index/star-fusion/v1.9/GRCh38_gencode_v33_CTAT_lib_Apr062020/ctat_genome_lib_build_dir/",
        "OPTIONS": "--CPU 12 --vis --examine_coding_effect"
    },
    "tabix_wes": {
        "OPTIONS": "-p vcf"
    },
    "vt_decompose": {
        "OPTIONS": "-s"
    },
    "vt_normalize": {
        "OPTIONS": ""
    },
    "wisecondor_predict": {
        "OPTIONS": "--plot --bed",
        "PON": "/mnt/beegfs/database/bioinfo/seqoia/original/annotations/Human/GRCh38/index/wisecondorx/reference_202003101247_5kb.npz"
    }
}

with open(snakemake.output["config"], "w") as config_stream:
    json_config = dumps(sample, indent=4)
    config_stream.write(json_config)
