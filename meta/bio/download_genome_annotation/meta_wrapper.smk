wildcard_constraints:
    gxf=r"gff3|gtf",


def local_lookup(dpath, within, default):
    """Return local(path) if path is within config"""
    tmp = lookup(dpath=dpath, within=within, default=None)
    if tmp is None:
        return default
    return local(tmp)


rule download_genome_annotation_reference_ensembl:
    """download genome gff"""
    output:
        temp(
            "temp/download_genome_annotation_reference_ensembl/{species}.{build}.{release}.{gxf}"
        ),
    log:
        "logs/download_genome_annotation_reference_ensembl/{species}.{build}.{release}.{gxf}.log",
    benchmark:
        "benchmarks/curl/download_genome_annotation_reference_ensembl_{species}.{build}.{release}.{gxf}.tsv"
    threads: 1
    params:
        species="{species}",
        build="{build}",
        release="{release}",
    wrapper:
        f"{config['repo']}/bio/reference/ensembl-annotation"


rule download_genome_annotation_agat_default_config:
    """define agat configuration for gff"""
    output:
        temp("temp/default_agat_config.yaml"),
    log:
        "logs/download_genome_annotation_agat_default_config.log",
    benchmark:
        "benchmarks/agat/download_genome_annotation_agat_default_config.tsv"
    group:
        "build_agat_config"
    shadow:
        "minimal"
    threads: 1
    params:
        command="config",
        extra="",
    wrapper:
        f"{config['repo']}/bio/agat"


rule download_genome_annotation_go_yq_update_gff_config:
    """define agat configuration for gff"""
    input:
        "temp/default_agat_config.yaml",
    output:
        "resources/{species}.{build}.{release}/annotations/{species}.{build}.{release}.agat_gff3_config.yaml",
    log:
        "logs/download_genome_annotation_go_yq_update_gff_config/{species}.{build}.{release}.gff3.log",
    benchmark:
        "benchmarks/go_yq/download_genome_annotation_go_yq_update_gff_config_{species}.{build}.{release}.gff3.tsv"
    group:
        "build_agat_config"
    threads: 1
    params:
        extra="--verbose",
        command="eval",
        expression=".cpu=6 | .log=false | .progress_bar=false | .minimum_chunk_size=1000",
    wrapper:
        f"{config['repo']}/utils/go-yq"


use rule download_genome_annotation_go_yq_update_gff_config as download_genome_annotation_go_yq_update_gtf_config with:
    """define agat configuration for gtf"""
    input:
        "resources/{species}.{build}.{release}/annotations/{species}.{build}.{release}.agat_gff3_config.yaml",
    output:
        "resources/{species}.{build}.{release}/annotations/{species}.{build}.{release}.agat_gtf_config.yaml",
    log:
        "logs/download_genome_annotation_go_yq_update_gtf_config/{species}.{build}.{release}.gtf.log",
    benchmark:
        "benchmarks/go_yq/download_genome_annotation_go_yq_update_gtf_config_{species}.{build}.{release}.gtf.tsv"
    group:
        "build_agat_config"
    threads: 1
    params:
        extra="--verbose",
        command="eval",
        expression='.output_format = "GTF"',


use rule download_genome_annotation_agat_default_config as download_genome_annotation_agat_convert_sp_gxf2gxf with:
    """fix common ensembl format issues"""
    input:
        gff="temp/download_genome_annotation_reference_ensembl/{species}.{build}.{release}.{gxf}",
        config="resources/{species}.{build}.{release}/annotations/{species}.{build}.{release}.agat_{gxf}_config.yaml",
    output:
        out=temp("temp/download_genome_annotation_agat_convert_sp_gxf2gxf/{species}.{build}.{release}.{gxf}"),
    log:
        "logs/download_genome_annotation_agat_convert_sp_gxf2gxf/{species}.{build}.{release}.{gxf}.log",
    benchmark:
        "benchmarks/agat/download_genome_annotation_agat_convert_sp_gxf2gxf_{species}.{build}.{release}.{gxf}.tsv"
    group:
        "agat_genome_annotation_processing"
    shadow:
        "minimal"
    threads: 6
    params:
        command="agat_convert_sp_gxf2gxf.pl",
        extra="",


use rule download_genome_annotation_agat_default_config as download_genome_annotation_agat_sq_filter_feature_from_fasta with:
    """ensure fasta and gff have the same contigs"""
    input:
        gff="temp/download_genome_annotation_agat_convert_sp_gxf2gxf/{species}.{build}.{release}.{gxf}",
        config="resources/{species}.{build}.{release}/annotations/{species}.{build}.{release}.agat_{gxf}_config.yaml",
        fasta=local_lookup(
            dpath="pathvars/genome_sequence",
            within=config,
            default="resources/{species}.{build}.{release}/sequences/{species}.{build}.{release}.dna.fasta",
        ),
    output:
        o=temp(
            "temp/download_genome_annotation_agat_sq_filter_feature_from_fasta/{species}.{build}.{release}.{gxf}"
        ),
    log:
        "logs/download_genome_annotation_agat_sq_filter_feature_from_fasta/{species}.{build}.{release}.{gxf}.log",
    benchmark:
        "benchmarks/agat/download_genome_annotation_agat_sq_filter_feature_from_fasta_{species}.{build}.{release}.{gxf}.tsv"
    group:
        "agat_genome_annotation_processing"
    shadow:
        "minimal"
    threads: 6
    params:
        extra="",
        command="agat_sq_filter_feature_from_fasta.pl",


use rule download_genome_annotation_agat_default_config as download_genome_annotation_agat_sp_filter_feature_by_attribute_value with:
    """filter out tsl na from genome annotation"""
    input:
        gff="temp/download_genome_annotation_agat_sq_filter_feature_from_fasta/{species}.{build}.{release}.{gxf}",
        config="resources/{species}.{build}.{release}/annotations/{species}.{build}.{release}.agat_{gxf}_config.yaml",
    output:
        o="resources/{species}.{build}.{release}/annotations/{species}.{build}.{release}.{gxf}",
    log:
        "logs/download_genome_annotation_agat_sp_filter_feature_by_attribute_value/{species}.{build}.{release}.{gxf}.log",
    benchmark:
        "benchmarks/agat/download_genome_annotation_agat_sp_filter_feature_by_attribute_value_{species}.{build}.{release}.{gxf}.tsv"
    group:
        "agat_genome_annotation_processing"
    shadow:
        "minimal"
    threads: 6
    params:
        extra="--attribute 'transcript_support_level' --value 'NA' --test '='",
        command="agat_sp_filter_feature_by_attribute_value.pl",


rule download_genome_annotation_sed_replace_genomes_data_in_report:
    """replace snakemake patterns for generic report"""
    input:
        workflow.source_path(
            "resources/reports/download_genome_annotation_material_and_methods.rst"
        ),
    output:
        temp(
            "temp/download_genome_annotation_sed_replace_genomes_data_in_report/material_and_methods_{species}_{build}_{release}.rst"
        ),
    log:
        "logs/download_genome_annotation_sed_replace_genomes_data_in_report/material_and_methods_{species}_{build}_{release}.log",
    benchmark:
        "benchmarks/sed/download_genome_annotation_sed_replace_genomes_data_in_report_{species}_{build}_{release}.tsv"
    group:
        "genome_annotation_material_and_methods"
    threads: 1
    params:
        extra="",
        expr="s| params.species |{species}|g;s| params.build |{build}|g;s| params.release |{release}|g;s|\\{\\{||g;s|\\}\\}||g",
    wrapper:
        f"{config['repo']}/utils/sed"


rule download_genome_annotation_material_and_methods:
    """provide comprehensible material and methods"""
    input:
        rst="temp/download_genome_annotation_sed_replace_genomes_data_in_report/material_and_methods_{species}_{build}_{release}.rst"
    output:
        html="resources/{species}.{build}.{release}/annotations/Genome_annotation_material_and_methods.html",
    log:
        "logs/download_genome_annotation_material_and_methods/{species}.{build}.{release}.log",
    benchmark:
        "benchmarks/rst2html/download_genome_annotation_material_and_methods_{species}.{build}.{release}.tsv"
    group:
        "genome_annotation_material_and_methods"
    threads: 1
    params:
        extra="--embed-stylesheet --title='Material and methods for genome annotation'",
    wrapper:
        f"{config['repo']}/utils/docutils/rst2html"
