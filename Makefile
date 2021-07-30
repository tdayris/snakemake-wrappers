SHELL := bash
.ONESHELL:
.SHELLFLAGS := -euic
.DELETE_ON_ERROR:
MAKEFLAGS += --warn-undefined-variables
MAKEFLAGS += --no-builtin-rules

### Variables ###
# Tools
PYTEST           = pytest
BASH             = bash
CONDA            = conda
PYTHON           = python3.8
SNAKEMAKE        = snakemake
CONDA_ACTIVATE   = source $$(conda info --base)/etc/profile.d/conda.sh && conda activate && conda activate


html:
	source $$(conda info --base)/etc/profile.d/conda.sh && \
	conda activate && \
	conda activate snakemake && \
	cd docs && \
	rm -rf wrappers/* && \
	rm -rf meta-wrappers/* && \
	rm -rf bigr_pipelines/* && \
	make clean && \
	make html && \
	cd - && \
	git add docs && \
	git commit -m "[doc] (HTML): Documentation update"
.PHONY: html


test:
	source $$(conda info --base)/etc/profile.d/conda.sh && \
	conda activate && \
	conda activate snakemake && \
	pytest -vvs 'test.py::$(ARGS)'
.PHONY: test

search:
	grep $(ARGS) test.py
