SHELL := bash
.ONESHELL:
.SHELLFLAGS := -eic
.DELETE_ON_ERROR:
MAKEFLAGS += --warn-undefined-variables
MAKEFLAGS += --no-builtin-rules

### Variables  ###
### Tools      ###

BASH             = bash
CONDA_ACTIVATE   = source "$$(conda info --base)/etc/profile.d/conda.sh" ; conda activate ; conda activate


### PATHS      ###

BASEDIR          = /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines

### Recipes    ###
### EaCoN      ###
cytoscan_eacon:
	cd ${BASEDIR}/cytoscan_eacon/test && \
	${BASH} ${BASEDIR}/cytoscan_eacon/run.sh && \
	rm -r config.yaml ctrl/ logs/ .snakemake/ && \
	cd -
.PHONY: cytoscan_eacon

oncoscan_eacon:
	cd ${BASEDIR}/oncoscan_eacon/test && \
	${BASH} ${BASEDIR}/oncoscan_eacon/run.sh && \
	rm -r config.yaml REF_20200513/ logs/ .snakemake/ && \
	cd -
.PHONY: oncoscan_eacon

### Salmon quantification ###
salmon_quant_hg19:
	cd ${BASEDIR}/salmon_quant/test && \
	${BASH} ${BASEDIR}/salmon_quant/run.sh hg19 && \
	rm -r fastp/ logs/ .snakemake/ multiqc/ reads/ salmon/ snakemake.salmon_quant.log config.yaml design.tsv && \
	cd -
.PHONY: salmon_quant_hg19

salmon_quant_hg38:
	cd ${BASEDIR}/salmon_quant/test && \
	${BASH} ${BASEDIR}/salmon_quant/run.sh hg38 && \
	rm -r fastp/ logs/ .snakemake/ multiqc/ reads/ salmon/ snakemake.salmon_quant.log config.yaml design.tsv && \
	cd -
.PHONY: salmon_quant_hg38

### DGE Deseq2 ###
dge_deseq2_hg19:
	cd ${BASEDIR}/dge_deseq2/test && \
	${BASH} ${BASEDIR}/dge_deseq2/run.sh hg19 && \
	rm -r config.yaml deseq2/ fastp/ figures/ logs/ multiqc/ reads/ .snakemake/ salmon/ tximport/ my_dds.RDS snakemake.dge_deseq2.log && \
	cd -
.PHONY: dge_deseq2_hg19

dge_deseq2_hg38:
	cd ${BASEDIR}/dge_deseq2/test && \
	${BASH} ${BASEDIR}/dge_deseq2/run.sh hg38 && \
	rm -r config.yaml deseq2/ fastp/ figures/ logs/ multiqc/ reads/ .snakemake/ salmon/ tximport/ my_dds.RDS snakemake.dge_deseq2.log && \
	cd -
.PHONY: dge_deseq2_hg38

### Ampliseq ###
variant_calling_ampliseq_hg19:
	cd ${BASEDIR}/variant_calling_ampliseq/test && \
	${BASH} ${BASEDIR}/variant_calling_ampliseq/run.sh hg19 && \
	rm -r bwa_mem2/ fastp/ gatk/ logs/ meta_caller/ multiqc/ mutect2/ reads/ snpeff/ snpsift/ varscan2/ snakemake.variant_calling_ampliseq.log .snakemake/ config_variant_calling_ampliseq.yaml config.yaml design.tsv && \
	cd -
.PHONY: variant_calling_ampliseq_hg19

variant_calling_ampliseq_hg38:
	cd ${BASEDIR}/variant_calling_ampliseq/test && \
	${BASH} ${BASEDIR}/variant_calling_ampliseq/run.sh hg38 && \
	rm -r bwa_mem2/ fastp/ gatk/ logs/ meta_caller/ multiqc/ mutect2/ reads/ snpeff/ snpsift/ varscan2/ snakemake.variant_calling_ampliseq.log .snakemake/ config_variant_calling_ampliseq.yaml config.yaml design.tsv && \
	cd -
.PHONY: variant_calling_ampliseq_hg38


### WES germline ###
variant_calling_germline_hg19:
	cd ${BASEDIR}/variant_calling_germline/test && \
	${BASH} ${BASEDIR}/variant_calling_germline/run.sh hg19 && \
	rm -r bwa_mem2/ fastp/ gatk/ logs/ meta_caller/ multiqc/ mutect2/ reads/ snpeff/ snpsift/ varscan2/ snakemake.variant_calling_germline.log .snakemake/ config_variant_calling_germline.yaml config.yaml design.tsv && \
	cd -
.PHONY: variant_calling_germline_hg19

variant_calling_germline_hg38:
	cd ${BASEDIR}/variant_calling_germline/test && \
	${BASH} ${BASEDIR}/variant_calling_germline/run.sh hg38 && \
	rm -r bwa_mem2/ fastp/ gatk/ logs/ meta_caller/ multiqc/ mutect2/ reads/ snpeff/ snpsift/ varscan2/ snakemake.variant_calling_germline.log .snakemake/ config_variant_calling_germline.yaml config.yaml design.tsv && \
	cd -
.PHONY: variant_calling_germline_hg38

### Fastqc_multiqc ###
fastqc_multiqc:
	cd ${BASEDIR}/fastqc_multiqc/test && \
	${BASH} ${BASEDIR}/fastqc_multiqc/run.sh && \
	cd -
.PHONY: fastqc_multiqc
