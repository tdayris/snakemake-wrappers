#!/usr/bin/python3.8
# -*- coding: utf-8 -*-

"""
This script creates an environment file used by
bigr_pipelines to guess and prepare design files,
config files, and launch the right pipeline with
the right parameters.

Command line example:
python3 pipeline_cache.py rnaseq --dge
"""

import argparse
import logging
import os
import sys
import yaml

from pathlib import Path
from snakemake.utils import makedirs

def parser() -> argparse.ArgumentParser:
    """
    Build the command line parser object
    """
    main_parser = argparse.ArgumentParser(
        description=sys.modules[__name__].__doc__,
        formatter_class=CustomFormatter,
        epilog="This script does any magic, please check results."
    )

    # Positional arguments
    main_parser.add_argument(
        "pipeline",
        help="Name of the pipeline you'd like to use",
        type=str,
        choices=["rnaseq", "wes_somatic", "none"],
    )

    # Optional arguments
    main_parser.add_argument(
        "--genome",
        type=str,
        help="Organism code",
        choices=["GRCh38", "GRCm38"],
    )

source_path = Path(os.path.realpath(__file__))  / ".." / ".." / ".."
exec_path = Path(".").absolute()
makedirs(".bigr_pipelines", mode=755, exist_ok=True)

cache_data = {
    "source": str(source_path.resolve()),
    "max_threads": 20,
    "design": str(exec_path / "data_input/"),
    "mm10": {
        "name": "GRCm38.99",
        "ncbi_build": "GRCm38",
        "patch": "99",
        "index": {
            "star": "/mnt/beegfs/database/bioinfo/Index_DB/STAR/GRCm38",
            "salmon": "/mnt/beegfs/database/bioinfo/Index_DB/Salmon/mm10/index",
            "msi": "/mnt/beegfs/database/bioinfo/Index_DB/MSI/GRCm38.99/GRCm38.99.msi.list",
            "bwa_index": [
                "/mnt/beegfs/database/bioinfo/Index_DB/BWA/2.2.1/GRCm38.99/GRCm38.99.mus_musculus.dna.fasta.0123",
                "/mnt/beegfs/database/bioinfo/Index_DB/BWA/2.2.1/GRCm38.99/GRCm38.99.mus_musculus.dna.fasta.amb",
                "/mnt/beegfs/database/bioinfo/Index_DB/BWA/2.2.1/GRCm38.99/GRCm38.99.mus_musculus.dna.fasta.ann",
                "/mnt/beegfs/database/bioinfo/Index_DB/BWA/2.2.1/GRCm38.99/GRCm38.99.mus_musculus.dna.fasta.bwt.2bit.64",
                "/mnt/beegfs/database/bioinfo/Index_DB/BWA/2.2.1/GRCm38.99/GRCm38.99.mus_musculus.dna.fasta.pac",
            ],
            "annot_sv": "/mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/annot_sv/source/AnnotSV/",
        },
        "fasta": {
            "genome": "/mnt/beegfs/database/bioinfo/Index_DB/Fasta/Ensembl/GRCm38.99/GRCm38.99.mus_musculus.dna.fasta",
            "genome_index": "/mnt/beegfs/database/bioinfo/Index_DB/Fasta/Ensembl/GRCm38.99/GRCm38.99.mus_musculus.dna.fasta.fai",
            "genome_dict": "/mnt/beegfs/database/bioinfo/Index_DB/Fasta/Ensembl/GRCm38.99/GRCm38.99.mus_musculus.dna.dict",
        },
        "annotations": {
            "gtf": "/mnt/beegfs/database/bioinfo/Index_DB/Fasta/Ensembl/GRCm38.99/GRCm38.99.mus_musculus.gtf",
            "reflat": "/mnt/beegfs/database/bioinfo/Index_DB/refFlat/GRCm38/refFlat.txt.gz",
            "refgene_model": "/mnt/beegfs/database/bioinfo/Index_DB/refgene/bed12/GRCm38/GRCm38.nochr.bed12",
            "af_only": None,
            "af_only_tbi": None,
            "dbsnp": "/mnt/beegfs/database/bioinfo/Index_DB/Fasta/Ensembl/GRCm38.99/GRCm38.99.mus_musculus.all.vcf.gz",
            "dbsnp_tbi": "/mnt/beegfs/database/bioinfo/Index_DB/Fasta/Ensembl/GRCm38.99/GRCm38.99.mus_musculus.all.vcf.gz.tbi",
        },
        "chromosomes": list(range(1, 23)) + list(map(str, range(1, 23))) + ["MT", "X", "Y"],
        "db": {
            "CTAT_resource_lib": "/mnt/beegfs/database/bioinfo/Index_DB/CTAT_Resource_Lib/GRCh38_gencode_v37/GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir/",
            "cibersort_mat": "/mnt/beegfs/software/cibersort/1.0.6/LM22.txt",
            "MSigDB_c1_Positional_Gene_Sets": "/mnt/beegfs/database/bioinfo/Index_DB/MSigDB/ENTREZID/7.5.1/c1.all.v7.5.1.entrez.gmt",
            "MSigDB_c2_Curated_Gene_Sets": "/mnt/beegfs/database/bioinfo/Index_DB/MSigDB/ENTREZID/7.5.1/c2.all.v7.5.1.entrez.gmt",
            "MSigDB_c3_Regulatory_Targets": "/mnt/beegfs/database/bioinfo/Index_DB/MSigDB/ENTREZID/7.5.1/c3.all.v7.5.1.entrez.gmt",
            "MSigDB_c4_Cancer_Oriented_MicroArray": "/mnt/beegfs/database/bioinfo/Index_DB/MSigDB/ENTREZID/7.5.1/c4.all.v7.5.1.entrez.gmt",
            "MSigDB_c5_Onthology_Gene_Sets": "/mnt/beegfs/database/bioinfo/Index_DB/MSigDB/ENTREZID/7.5.1/c5.all.v7.5.1.entrez.gmt",
            "MSigDB_c6_Oncogenic_Signatures": "/mnt/beegfs/database/bioinfo/Index_DB/MSigDB/ENTREZID/7.5.1/c6.all.v7.5.1.entrez.gmt",
            "MSigDB_c7_Immunologic_Signatures": "/mnt/beegfs/database/bioinfo/Index_DB/MSigDB/ENTREZID/7.5.1/c7.all.v7.5.1.entrez.gmt",
            "MSigDB_c8_Cell_Type_Signatures": "/mnt/beegfs/database/bioinfo/Index_DB/MSigDB/ENTREZID/7.5.1/c8.all.v7.5.1.entrez.gmt",
            "MSigDB_hallmark": "/mnt/beegfs/database/bioinfo/Index_DB/MSigDB/ENTREZID/7.5.1/h.all.v7.5.1.entrez.gmt",
            "Chemical_genetic_perturbation": "/mnt/beegfs/database/bioinfo/Index_DB/MSigDB/ENTREZID/7.5.1/c2.cgp.v7.5.1.entrez.gmt",
            "Canonical_pathways": "/mnt/beegfs/database/bioinfo/Index_DB/MSigDB/ENTREZID/7.5.1/c2.cp.v7.5.1.entrez.gmt",
            "Biocarta": "/mnt/beegfs/database/bioinfo/Index_DB/MSigDB/ENTREZID/7.5.1/c2.cp.biocarta.v7.5.1.entrez.gmt",
            "KEGG": "/mnt/beegfs/database/bioinfo/Index_DB/MSigDB/ENTREZID/7.5.1/c2.cp.kegg.v7.5.1.entrez.gmt",
            "PID": "/mnt/beegfs/database/bioinfo/Index_DB/MSigDB/ENTREZID/7.5.1/c2.cp.pid.v7.5.1.entrez.gmt",
            "MSigDB_Reactome": "/mnt/beegfs/database/bioinfo/Index_DB/MSigDB/ENTREZID/7.5.1/c2.cp.reactome.v7.5.1.entrez.gmt",
            "MSigDB_WikiPathways": "/mnt/beegfs/database/bioinfo/Index_DB/MSigDB/ENTREZID/7.5.1/c2.cp.wikipathways.v7.5.1.entrez.gmt",
            "MicroRNA_targets": "/mnt/beegfs/database/bioinfo/Index_DB/MSigDB/ENTREZID/7.5.1/c3.mir.v7.5.1.entrez.gmt",
            "miRDB": "/mnt/beegfs/database/bioinfo/Index_DB/MSigDB/ENTREZID/7.5.1/c3.mir.mirdb.v7.5.1.entrez.gmt",
            "MicroRNA_legacy": "/mnt/beegfs/database/bioinfo/Index_DB/MSigDB/ENTREZID/7.5.1/c3.mir.mir_legacy.v7.5.1.entrez.gmt",
            "Transcription_factor_targets": "/mnt/beegfs/database/bioinfo/Index_DB/MSigDB/ENTREZID/7.5.1/c3.tft.v7.5.1.entrez.gmt",
            "Transcription_factor_GTRD_Kolmykov": "/mnt/beegfs/database/bioinfo/Index_DB/MSigDB/ENTREZID/7.5.1/c3.tft.gtrd.v7.5.1.entrez.gmt",
            "Transcription_factor_targets_legacy": "/mnt/beegfs/database/bioinfo/Index_DB/MSigDB/ENTREZID/7.5.1/c3.tft.tft_legacy.v7.5.1.entrez.gmt",
            "Cancer_gene_neighborhood": "/mnt/beegfs/database/bioinfo/Index_DB/MSigDB/ENTREZID/7.5.1/c4.cgn.v7.5.1.entrez.gmt",
            "Cancer_modules": "/mnt/beegfs/database/bioinfo/Index_DB/MSigDB/ENTREZID/7.5.1/c4.cm.v7.5.1.entrez.gmt",
            "Gene_Othology_all": "/mnt/beegfs/database/bioinfo/Index_DB/MSigDB/ENTREZID/7.5.1/c5.go.v7.5.1.entrez.gmt",
            "Gene_Othology_BP": "/mnt/beegfs/database/bioinfo/Index_DB/MSigDB/ENTREZID/7.5.1/c5.go.bp.v7.5.1.entrez.gmt",
            "Gene_Othology_CC": "/mnt/beegfs/database/bioinfo/Index_DB/MSigDB/ENTREZID/7.5.1/c5.go.cc.v7.5.1.entrez.gmt",
            "Gene_Othology_MF": "/mnt/beegfs/database/bioinfo/Index_DB/MSigDB/ENTREZID/7.5.1/c5.go.mf.v7.5.1.entrez.gmt",
            "MSigDB_Human_Phenotype_Onthology": "/mnt/beegfs/database/bioinfo/Index_DB/MSigDB/ENTREZID/7.5.1/c5.hpo.v7.5.1.entrez.gmt",
            "ImmuneSigDB": "/mnt/beegfs/database/bioinfo/Index_DB/MSigDB/ENTREZID/7.5.1/c7.immunesigdb.v7.5.1.entrez.gmt",
            "Vaccine_Response": "/mnt/beegfs/database/bioinfo/Index_DB/MSigDB/ENTREZID/7.5.1/c7.vax.v7.5.1.entrez.gmt",
            "WikiPathway": "/mnt/beegfs/database/bioinfo/Index_DB/WikiPathway/HomoSapiens/hsapiens.WP.ENSG.gmt",
            "Reactome": "/mnt/beegfs/database/bioinfo/Index_DB/Reactome/HomoSapiens/hsapiens.REAC.ENSG.gmt",
            "Human_Phenotype_Onthology": "/mnt/beegfs/database/bioinfo/Index_DB/HumanPhenotypeOnthology/HomoSapiens/hsapiens.HP.ENSG.gmt",
            "Human_Protein_Atlas": "/mnt/beegfs/database/bioinfo/Index_DB/HumanProteinAtlas/HomoSapiens/hsapiens.HPA.ENSG.gmt",
            "MiRNABase": "/mnt/beegfs/database/bioinfo/Index_DB/MiRNABase/HomoSapiens/hsapiens.MIRNA.ENSG.gmt",
            "CORUM": "/mnt/beegfs/database/bioinfo/Index_DB/CORUM/HomoSapiens/hsapiens.CORUM.ENSG.gmt",
            "SMART": "/mnt/beegfs/database/bioinfo/Index_DB/PPI_Hotspot/PPI_Protein_Domains__SMART_.tsv",
            "Pfam": "/mnt/beegfs/database/bioinfo/Index_DB/PPI_Hotspot/PPI_Protein_Domains__Pfam_.tsv",
            "InterPro": "/mnt/beegfs/database/bioinfo/Index_DB/PPI_Hotspot/PPI_Protein_Domains_and_Features__InterPro_.tsv",
            "Tissues": "/mnt/beegfs/database/bioinfo/Index_DB/PPI_Hotspot/PPI_Tissue_expression__TISSUES_.tsv",
            "STRING": "/mnt/beegfs/database/bioinfo/Index_DB/PPI_Hotspot/PPI_Local_Network_Cluster__STRING_.tsv",
            "Disease": "/mnt/beegfs/database/bioinfo/Index_DB/PPI_Hotspot/PPI_Disease-gene_associations__DISEASES_.tsv",
            "snpeff": "/mnt/beegfs/database/bioinfo/Index_DB/SnpEff/GRCm38.86",
            "kaviar": None,
            "kaviar_tbi": None,
            "cosmic": None,
            "cosmic_tbi": None,
            "dbnsfp": None,
            "clinvar": None,
            "clinvar_tbi": None,
            "gnomad": None,
            "gnomad_tbi": None,
            "exac": None,
            "exac_tbi": None,
            "dbvar": None,
            "dbvar_tbi": None,
            "gwascat": None,
            "mane": None,
            "revel": None,
            "mistic": None,
            "cancer_census": None,
            "oncokb": None,
        }
    },
    "hg38": {
        "name": "GRCh38.99",
        "ncbi_build": "GRCh38",
        "patch": "99",
        "index": {
            "star": "/mnt/beegfs/database/bioinfo/Index_DB/STAR/GRCh38",
            "salmon": "/mnt/beegfs/database/bioinfo/Index_DB/Salmon/hg38/index",
            "msi": "/mnt/beegfs/database/bioinfo/Index_DB/MSI/GRCh38.99/GRCh38.99.msi.list",
            "bwa_index": [
                "/mnt/beegfs/database/bioinfo/Index_DB/BWA/2.2.1/GRCh38.99/GRCh38.99.homo_sapiens.dna.fasta.0123",
                "/mnt/beegfs/database/bioinfo/Index_DB/BWA/2.2.1/GRCh38.99/GRCh38.99.homo_sapiens.dna.fasta.amb",
                "/mnt/beegfs/database/bioinfo/Index_DB/BWA/2.2.1/GRCh38.99/GRCh38.99.homo_sapiens.dna.fasta.ann",
                "/mnt/beegfs/database/bioinfo/Index_DB/BWA/2.2.1/GRCh38.99/GRCh38.99.homo_sapiens.dna.fasta.bwt.2bit.64",
                "/mnt/beegfs/database/bioinfo/Index_DB/BWA/2.2.1/GRCh38.99/GRCh38.99.homo_sapiens.dna.fasta.pac",
            ],
            "annot_sv": "/mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/annot_sv/source/AnnotSV/",
        },
        "fasta": {
            "genome": "/mnt/beegfs/database/bioinfo/Index_DB/Fasta/Ensembl/GRCh38.99/GRCh38.99.homo_sapiens.dna.main_chr.fasta",
            "genome_index": "/mnt/beegfs/database/bioinfo/Index_DB/Fasta/Ensembl/GRCh38.99/GRCh38.99.homo_sapiens.dna.main_chr.fasta.fai",
            "genome_dict": "/mnt/beegfs/database/bioinfo/Index_DB/Fasta/Ensembl/GRCh38.99/GRCh38.99.homo_sapiens.dna.main_chr.dict"
        },
        "annotations": {
            "gtf": "/mnt/beegfs/database/bioinfo/Index_DB/Fasta/Ensembl/GRCh38.99/GRCh38.99.homo_sapiens.gtf",
            "reflat": "/mnt/beegfs/database/bioinfo/Index_DB/refFlat/GRCh38/refFlat.txt.gz",
            "refgene_model": "/mnt/beegfs/database/bioinfo/Index_DB/refgene/bed12/GRCh38/GRCh38.nochr.bed12",
            "af_only":"/mnt/beegfs/database/bioinfo/Index_DB/GATK/mutect2_gnomad_af_only/hg38/somatic-hg38_af-only-gnomad.hg38.nochr.vcf.gz",
            "af_only_tbi": "/mnt/beegfs/database/bioinfo/Index_DB/GATK/mutect2_gnomad_af_only/hg38/somatic-hg38_af-only-gnomad.hg38.nochr.vcf.gz.tbi",
            "dbsnp": "/mnt/beegfs/database/bioinfo/Index_DB/dbSNP/GRCh38p7/All_20180418.nochr.vcf.gz",
            "dbsnp_tbi": "/mnt/beegfs/database/bioinfo/Index_DB/dbSNP/GRCh38p7/All_20180418.nochr.vcf.gz.tbi",
        },
        "chromosomes": list(range(1, 23)) + list(map(str, range(1, 23))) + ["MT", "X", "Y"],
        "db": {
            "CTAT_resource_lib": "/mnt/beegfs/database/bioinfo/Index_DB/CTAT_Resource_Lib/GRCh38_gencode_v37/GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir/",
            "cibersort_mat": "/mnt/beegfs/software/cibersort/1.0.6/LM22.txt",
            "MSigDB_c1_Positional_Gene_Sets": "/mnt/beegfs/database/bioinfo/Index_DB/MSigDB/ENTREZID/7.5.1/c1.all.v7.5.1.entrez.gmt",
            "MSigDB_c2_Curated_Gene_Sets": "/mnt/beegfs/database/bioinfo/Index_DB/MSigDB/ENTREZID/7.5.1/c2.all.v7.5.1.entrez.gmt",
            "MSigDB_c3_Regulatory_Targets": "/mnt/beegfs/database/bioinfo/Index_DB/MSigDB/ENTREZID/7.5.1/c3.all.v7.5.1.entrez.gmt",
            "MSigDB_c4_Cancer_Oriented_MicroArray": "/mnt/beegfs/database/bioinfo/Index_DB/MSigDB/ENTREZID/7.5.1/c4.all.v7.5.1.entrez.gmt",
            "MSigDB_c5_Onthology_Gene_Sets": "/mnt/beegfs/database/bioinfo/Index_DB/MSigDB/ENTREZID/7.5.1/c5.all.v7.5.1.entrez.gmt",
            "MSigDB_c6_Oncogenic_Signatures": "/mnt/beegfs/database/bioinfo/Index_DB/MSigDB/ENTREZID/7.5.1/c6.all.v7.5.1.entrez.gmt",
            "MSigDB_c7_Immunologic_Signatures": "/mnt/beegfs/database/bioinfo/Index_DB/MSigDB/ENTREZID/7.5.1/c7.all.v7.5.1.entrez.gmt",
            "MSigDB_c8_Cell_Type_Signatures": "/mnt/beegfs/database/bioinfo/Index_DB/MSigDB/ENTREZID/7.5.1/c8.all.v7.5.1.entrez.gmt",
            "MSigDB_hallmark": "/mnt/beegfs/database/bioinfo/Index_DB/MSigDB/ENTREZID/7.5.1/h.all.v7.5.1.entrez.gmt",
            "Chemical_genetic_perturbation": "/mnt/beegfs/database/bioinfo/Index_DB/MSigDB/ENTREZID/7.5.1/c2.cgp.v7.5.1.entrez.gmt",
            "Canonical_pathways": "/mnt/beegfs/database/bioinfo/Index_DB/MSigDB/ENTREZID/7.5.1/c2.cp.v7.5.1.entrez.gmt",
            "Biocarta": "/mnt/beegfs/database/bioinfo/Index_DB/MSigDB/ENTREZID/7.5.1/c2.cp.biocarta.v7.5.1.entrez.gmt",
            "KEGG": "/mnt/beegfs/database/bioinfo/Index_DB/MSigDB/ENTREZID/7.5.1/c2.cp.kegg.v7.5.1.entrez.gmt",
            "PID": "/mnt/beegfs/database/bioinfo/Index_DB/MSigDB/ENTREZID/7.5.1/c2.cp.pid.v7.5.1.entrez.gmt",
            "MSigDB_Reactome": "/mnt/beegfs/database/bioinfo/Index_DB/MSigDB/ENTREZID/7.5.1/c2.cp.reactome.v7.5.1.entrez.gmt",
            "MSigDB_WikiPathways": "/mnt/beegfs/database/bioinfo/Index_DB/MSigDB/ENTREZID/7.5.1/c2.cp.wikipathways.v7.5.1.entrez.gmt",
            "MicroRNA_targets": "/mnt/beegfs/database/bioinfo/Index_DB/MSigDB/ENTREZID/7.5.1/c3.mir.v7.5.1.entrez.gmt",
            "miRDB": "/mnt/beegfs/database/bioinfo/Index_DB/MSigDB/ENTREZID/7.5.1/c3.mir.mirdb.v7.5.1.entrez.gmt",
            "MicroRNA_legacy": "/mnt/beegfs/database/bioinfo/Index_DB/MSigDB/ENTREZID/7.5.1/c3.mir.mir_legacy.v7.5.1.entrez.gmt",
            "Transcription_factor_targets": "/mnt/beegfs/database/bioinfo/Index_DB/MSigDB/ENTREZID/7.5.1/c3.tft.v7.5.1.entrez.gmt",
            "Transcription_factor_GTRD_Kolmykov": "/mnt/beegfs/database/bioinfo/Index_DB/MSigDB/ENTREZID/7.5.1/c3.tft.gtrd.v7.5.1.entrez.gmt",
            "Transcription_factor_targets_legacy": "/mnt/beegfs/database/bioinfo/Index_DB/MSigDB/ENTREZID/7.5.1/c3.tft.tft_legacy.v7.5.1.entrez.gmt",
            "Cancer_gene_neighborhood": "/mnt/beegfs/database/bioinfo/Index_DB/MSigDB/ENTREZID/7.5.1/c4.cgn.v7.5.1.entrez.gmt",
            "Cancer_modules": "/mnt/beegfs/database/bioinfo/Index_DB/MSigDB/ENTREZID/7.5.1/c4.cm.v7.5.1.entrez.gmt",
            "Gene_Othology_all": "/mnt/beegfs/database/bioinfo/Index_DB/MSigDB/ENTREZID/7.5.1/c5.go.v7.5.1.entrez.gmt",
            "Gene_Othology_BP": "/mnt/beegfs/database/bioinfo/Index_DB/MSigDB/ENTREZID/7.5.1/c5.go.bp.v7.5.1.entrez.gmt",
            "Gene_Othology_CC": "/mnt/beegfs/database/bioinfo/Index_DB/MSigDB/ENTREZID/7.5.1/c5.go.cc.v7.5.1.entrez.gmt",
            "Gene_Othology_MF": "/mnt/beegfs/database/bioinfo/Index_DB/MSigDB/ENTREZID/7.5.1/c5.go.mf.v7.5.1.entrez.gmt",
            "MSigDB_Human_Phenotype_Onthology": "/mnt/beegfs/database/bioinfo/Index_DB/MSigDB/ENTREZID/7.5.1/c5.hpo.v7.5.1.entrez.gmt",
            "ImmuneSigDB": "/mnt/beegfs/database/bioinfo/Index_DB/MSigDB/ENTREZID/7.5.1/c7.immunesigdb.v7.5.1.entrez.gmt",
            "Vaccine_Response": "/mnt/beegfs/database/bioinfo/Index_DB/MSigDB/ENTREZID/7.5.1/c7.vax.v7.5.1.entrez.gmt",
            "WikiPathway": "/mnt/beegfs/database/bioinfo/Index_DB/WikiPathway/HomoSapiens/hsapiens.WP.ENSG.gmt",
            "Reactome": "/mnt/beegfs/database/bioinfo/Index_DB/Reactome/HomoSapiens/hsapiens.REAC.ENSG.gmt",
            "Human_Phenotype_Onthology": "/mnt/beegfs/database/bioinfo/Index_DB/HumanPhenotypeOnthology/HomoSapiens/hsapiens.HP.ENSG.gmt",
            "Human_Protein_Atlas": "/mnt/beegfs/database/bioinfo/Index_DB/HumanProteinAtlas/HomoSapiens/hsapiens.HPA.ENSG.gmt",
            "MiRNABase": "/mnt/beegfs/database/bioinfo/Index_DB/MiRNABase/HomoSapiens/hsapiens.MIRNA.ENSG.gmt",
            "CORUM": "/mnt/beegfs/database/bioinfo/Index_DB/CORUM/HomoSapiens/hsapiens.CORUM.ENSG.gmt",
            "SMART": "/mnt/beegfs/database/bioinfo/Index_DB/PPI_Hotspot/PPI_Protein_Domains__SMART_.tsv",
            "Pfam": "/mnt/beegfs/database/bioinfo/Index_DB/PPI_Hotspot/PPI_Protein_Domains__Pfam_.tsv",
            "InterPro": "/mnt/beegfs/database/bioinfo/Index_DB/PPI_Hotspot/PPI_Protein_Domains_and_Features__InterPro_.tsv",
            "Tissues": "/mnt/beegfs/database/bioinfo/Index_DB/PPI_Hotspot/PPI_Tissue_expression__TISSUES_.tsv",
            "STRING": "/mnt/beegfs/database/bioinfo/Index_DB/PPI_Hotspot/PPI_Local_Network_Cluster__STRING_.tsv",
            "Disease": "/mnt/beegfs/database/bioinfo/Index_DB/PPI_Hotspot/PPI_Disease-gene_associations__DISEASES_.tsv",
            "snpeff": "/mnt/beegfs/database/bioinfo/Index_DB/SnpEff/GRCh38.86",
            "kaviar": "/mnt/beegfs/database/bioinfo/Index_DB/Kaviar/HG38/Kaviar-160204-Public/vcfs/Kaviar-160204-Public-hg38-trim.vcf.gz",
            "kaviar_tbi": "/mnt/beegfs/database/bioinfo/Index_DB/Kaviar/HG38/Kaviar-160204-Public/vcfs/Kaviar-160204-Public-hg38-trim.vcf.gz.tbi",
            "cosmic": "/mnt/beegfs/database/bioinfo/COSMIC/73_20150629/CosmicCodingMuts.ucsc.nochr.vcf.gz",
            "cosmic_tbi": "/mnt/beegfs/database/bioinfo/COSMIC/73_20150629/CosmicCodingMuts.ucsc.nochr.vcf.gz.tbi",
            "dbnsfp": "/mnt/beegfs/database/bioinfo/Index_DB/dbNSFP/4.1/GRCh38/dbNSFP4.1a.txt.gz",
            "clinvar": "/mnt/beegfs/database/bioinfo/Index_DB/ClinVar/GRCh38/clinvar_20210404.vcf.gz",
            "clinvar_tbi": "/mnt/beegfs/database/bioinfo/Index_DB/ClinVar/GRCh38/clinvar_20210404.vcf.gz.tbi",
            "gnomad": "/mnt/beegfs/database/bioinfo/Index_DB/gnomad/gnomad.genomes.r3.0.sites.vcf.nochr.bgz",
            "gnomad_tbi": "/mnt/beegfs/database/bioinfo/Index_DB/gnomad/gnomad.genomes.r3.0.sites.vcf.nochr.bgz.tbi",
            "exac": "/mnt/beegfs/database/bioinfo/Index_DB/Exac/release1/ExAC.r1.sites.vep.fixed.vcf.gz",
            "exac_tbi": "/mnt/beegfs/database/bioinfo/Index_DB/Exac/release1/ExAC.r1.sites.vep.fixed.vcf.gz.tbi",
            "dbvar": "/mnt/beegfs/database/bioinfo/Index_DB/dbVar/GRCh38.variant_call.all.vcf.gz",
            "dbvar_tbi": "/mnt/beegfs/database/bioinfo/Index_DB/dbVar/GRCh38.variant_call.all.vcf.gz.tbi",
            "gwascat": "/mnt/beegfs/database/bioinfo/Index_DB/GWASCatalog/gwas_catalog_v1.0.2-studies_r2020-05-03.tsv",
            "mane": "/mnt/beegfs/database/bioinfo/Index_DB/MANE/MANE_human/release_0.95/MANE.GRCh38.v0.95.summary.txt.gz",
            "revel": "/mnt/beegfs/database/bioinfo/Index_DB/REVEL/revel_with_transcript_ids",
            "mistic": "/mnt/beegfs/database/bioinfo/Index_DB/MISTIC/MISTIC_GRCh38.tsv.gz",
            "cancer_census": "/mnt/beegfs/database/bioinfo/Index_DB/CancerGeneCensus/Census_allTue_Aug_31_15_11_39_2021.csv",
            "oncokb": "/mnt/beegfs/database/bioinfo/Index_DB/OncoKB/OncoKB.csv",
        }
    }
}