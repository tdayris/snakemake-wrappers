import logging
import os
import pandas
import sys
from pathlib import Path

logging.basicConfig(
    filename="snakemake.sigprofiler.log",
    filemode="w",
    level=logging.DEBUG
)

worflow_source_dir = Path(next(iter(workflow.get_sources()))).absolute().parent
common = str(worflow_source_dir / "../common/python")
sys.path.append(common)

from file_manager import *
from files_linker import *
from write_yaml import *
from messages import *
from snakemake.utils import min_version
min_version("6.0")

default_config = read_yaml(worflow_source_dir / "config.hg38.yaml")
configfile: get_config(default_config)
design = get_design(os.getcwd(), search_vcf_files).head(1)
design["Sample_id"] = design["Sample_id"].str.replace("-", "_")

install_genome = False

container: "docker://continuumio/miniconda3:4.4.10"
#localrules: bigr_copy

samples_list = design["Sample_id"].tolist()

wildcard_constraints:
    sample = r"|".join(samples_list)

vcf_links = link_vcf(
    design.Sample_id,
    design.Upstream_file
)

organism = (config.get("params", {"genome_build": "GRCh38"})
                  .get("genome_build", "GRCh38"))

rule target:
    input:
        expand(
            "SigProfiler/{sample}/test/input/",
            "SigProfiler/{sample}/Res/JOB_METADATA.txt",
            "SigProfiler/{sample}/test/output",
            #"DBS/{sample}/test/output",
            #"ID/{sample}/test/output",
            sample=samples_list
        )


rule sigprofiler_single_sample_sbs:
    input:
        vcf = "data_input/calls/{sample}.vcf.gz",
        sigprofiler_script = "/mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/sigprofiler_signatures/scripts/single_sample.py"
    output:
        outdir = directory("SigProfiler/{sample}/test/output"),
        vcf = "SigProfiler/{sample}/test/{sample}.vcf"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 4,
        time_min=lambda wildcards, attempt: attempt * 45,
        tmpdir="tmp"
    log:
        "logs/single_sample/{sample}.log"
    conda:
        "/mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/sigprofiler_signatures/env/single_sample_sigprofiler.yaml"
    params:
        mk = lambda wildcards: f"--parents --verbose SigProfiler/{wildcards.sample}/test/",
        gz = "--stdout --force",
        org = organism,
        gr = '-vP "^#"',
        cut = "-f 1-5",
        sed = "'/^chr/! s/^/chr/g'",
        install = "--install" if install_genome is True else ""
    shell:
        "mkdir {params.mk} > {log} 2>&1 && "
        "(gunzip {params.gz} {input.vcf} | "
        " grep {params.gr} | "
        " cut {params.cut} | "
        " sed {params.sed} ) > {output.vcf} 2>> {log} && "
        "python3 {input.sigprofiler_script} {output.vcf} {output.outdir}"
        " --organism {params.org} {params.install} >> {log} 2>&1 "



rule sigprofiler_single_sample_dbs:
    input:
        vcf = "data_input/calls/{sample}.vcf.gz",
        sigprofiler_script = "/mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/sigprofiler_signatures/scripts/single_sample.py",
        dbs = "/mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/conda/e31318e4e678218a151c58b51adefff7/lib/python3.7/site-packages/sigproSS/input/sigProfiler_DBS_signatures.csv"
    output:
        outdir = directory("DBS/{sample}/test/output"),
        vcf = "DBS/{sample}/test/{sample}.vcf"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 4,
        time_min=lambda wildcards, attempt: attempt * 45,
        tmpdir="tmp"
    log:
        "logs/single_sample/{sample}.log"
    conda:
        "/mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/sigprofiler_signatures/env/single_sample_sigprofiler.yaml"
    params:
        mk = lambda wildcards: f"--parents --verbose DBS/{wildcards.sample}/test/",
        gz = "--stdout --force",
        org = organism,
        gr = '-vP "^#"',
        cut = "-f 1-5",
        sed = "'/^chr/! s/^/chr/g'",
        install = "--install" if install_genome is True else "",
    shell:
        "mkdir {params.mk} > {log} 2>&1 && "
        "(gunzip {params.gz} {input.vcf} | "
        " grep {params.gr} | "
        " cut {params.cut} | "
        " sed {params.sed} ) > {output.vcf} 2>> {log} && "
        "python3 {input.sigprofiler_script} {output.vcf} {output.outdir}"
        " --dbs {input.dbs}"
        " --organism {params.org} {params.install} >> {log} 2>&1 "


rule sigprofiler_single_sample_id:
    input:
        vcf = "data_input/calls/{sample}.vcf.gz",
        sigprofiler_script = "/mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/sigprofiler_signatures/scripts/single_sample.py",
        id = "/mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/conda/e31318e4e678218a151c58b51adefff7/lib/python3.7/site-packages/sigproSS/input/sigProfiler_ID_signatures.csv"
    output:
        outdir = directory("ID/{sample}/test/output"),
        vcf = "ID/{sample}/test/{sample}.vcf"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 4,
        time_min=lambda wildcards, attempt: attempt * 45,
        tmpdir="tmp"
    log:
        "logs/single_sample/{sample}.log"
    conda:
        "/mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/sigprofiler_signatures/env/single_sample_sigprofiler.yaml"
    params:
        mk = lambda wildcards: f"--parents --verbose ID/{wildcards.sample}/test/",
        gz = "--stdout --force",
        org = organism,
        gr = '-vP "^#"',
        cut = "-f 1-5",
        sed = "'/^chr/! s/^/chr/g'",
        install = "--install" if install_genome is True else "",
    shell:
        "mkdir {params.mk} > {log} 2>&1 && "
        "(gunzip {params.gz} {input.vcf} | "
        " grep {params.gr} | "
        " cut {params.cut} | "
        " sed {params.sed} ) > {output.vcf} 2>> {log} && "
        "python3 {input.sigprofiler_script} {output.vcf} {output.outdir}"
        " --id {input.id}"
        " --organism {params.org} {params.install} >> {log} 2>&1 "


rule sigprofiler_matrix_generator:
    input:
        vcf = "data_input/calls/{sample}.vcf.gz",
        sigprofiler_script = "/mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/sigprofiler_signatures/scripts/matrix_generator.py"
    output:
        vcf = "SigProfiler/{sample}/test/{sample}.vcf",
        data_input = directory("SigProfiler/{sample}/test/input/"),
        logs = directory("SigProfiler/{sample}/test/logs"),
        out = directory(expand("SigProfiler/{sample}/test/output/{signatures}", signatures=["DBS", "ID", "SBS", "plots"], sample="{sample}")),
        out_vcf = directory(expand("SigProfiler/{sample}/test/output/vcf_files/{signatures}", signatures=["DBS", "ID", "SNV"], sample="{sample}"))
    message:
        "Building substitution matrices with SigProfiler on {wildcards.sample}"
    threads: 4
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 3,
        time_min=lambda wildcards, attempt: attempt * 60,
        tmpdir="tmp"
    conda:
        "/mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/sigprofiler_signatures/env/matrix_generator.yaml"
    log:
        "logs/sigprofiler/matrix_generator/{sample}.log"
    params:
        mk = lambda wildcards: f"--parents --verbose SigProfiler/{wildcards.sample}/test/",
        gz = "--stdout --force",
        org = organism,
        gr = '-vP "^#"',
        cut = "-f 1-5",
        sed = "'/^chr/! s/^/chr/g'",
        install = "--install" if install_genome is True else ""
    shell:
        "mkdir {params.mk} > {log} 2>&1 && "
        "(gunzip {params.gz} {input.vcf} | "
        " grep {params.gr} | "
        " cut {params.cut} | "
        " sed {params.sed} ) > {output.vcf} 2>> {log} && "
        "python3 {input.sigprofiler_script} {output.vcf} "
        " --organism {params.org} --sample-name {wildcards.sample} {params.install} "
        ">> {log} 2>&1 "


rule signature_extractor:
    input:
        matrices = "SigProfiler/{sample}/test/input",
        logs = "SigProfiler/{sample}/test/logs",
        out = expand("SigProfiler/{sample}/test/output/{signatures}", signatures=["DBS", "ID", "SBS", "plots"], sample="{sample}"),
        out_vcf = expand("SigProfiler/{sample}/test/output/vcf_files/{signatures}", signatures=["DBS", "ID", "SNV"], sample="{sample}"),
        sigprofiler_script = "/mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/sigprofiler_signatures/scripts/signature_extractor.py"
    output:
        sbs = directory("SigProfiler/{sample}/test/Res/SBS96/"),
        dsb = directory("SigProfiler/{sample}/test/Res/DBS78/"),
        job = "SigProfiler/{sample}/Res/JOB_METADATA.txt",
        seeds = "SigProfiler/{sample}/Res/Seeds.txt"
    message:
        "Extracting signatures in {wildcards.sample}"
    threads: 10
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 7,
        time_min=lambda wildcards, attempt: attempt * 60,
        tmpdir="tmp",
        #gres="gpu:t4:1"
    log:
        "logs/sigprofiler/extractor/{sample}.log"
    conda:
        "/mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/sigprofiler_signatures/env/extractor.yaml"
    params:
        sig = lambda wildcards: f"SigProfiler/{wildcards.sample}",
        org = organism,
        install = "--install" if install_genome is True else ""
    shell:
        "python3 {input.sigprofiler_script} {params.sig}/test {params.sig}/Res --threads {threads} --organism {params.org} {params.install} > {log} 2>&1"


rule sigprofiler_decompose:
    input:
        sbs = "SigProfiler/{sample}/Res/SBS96/",
        dsb = "SigProfiler/{sample}/Res/DBS78/",
        sigprofiler_script = "/mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/sigprofiler_signatures/scripts/sigprofiler_decompose.py"
    output:
        directory("SigProfiler/{sample}/Res/SBS96/Deconvolution_SB96_DeNovo")
    message:
        "Decomposing sugnatures on {wildcards.sample}"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
        time_min=lambda wildcards, attempt: attempt * 15,
        tmpdir="tmp"
    log:
        "logs/sigprofiler/decompose/{sample}.log"
    params:
        organism = organism,
        rs = "--verbose --checksum --recursive --human-readable --update",
        mk = "--parents --verbose",
        install = "--install" if install_genome is True else ""
    conda:
        "/mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/sigprofiler_signatures/env/extractor.yaml"
    shell:
        "rsync {params.rs} SigProfiler/{wildcards.sample}/test/ {output} > {log} 2>&1 && "
        "mkdir {params.mk} {output}/Res/SBS96/Suggested_Solution/Decomposed_Solution/ >> {log} 2>&1 && "
        "mkdir {params.mk} {output}/Res/SBS96/Suggested_Solution/De_Novo_Solution/ >> {log} 2>&1 && "
        "mkdir {params.mk} {output}/Res/ID83/Suggested_Solution/Decomposed_Solution/ >> {log} 2>&1 && "
        "mkdir {params.mk} {output}/Res/ID83/Suggested_Solution/De_Novo_Solution/ >> {log} 2>&1 && "
        "mkdir {params.mk} {output}/Res/DBS78/Suggested_Solution/Decomposed_Solution/ >> {log} 2>&1 && "
        "mkdir {params.mk} {output}/Res/DBS78/Suggested_Solution/De_Novo_Solution/ >> {log} 2>&1 && "
        "python3 {input.sigprofiler_script} {output} --organism {params.organism} --verbose {params.install} >> {log} 2>&1"