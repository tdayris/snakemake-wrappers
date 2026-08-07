wildcard_constraints:
    gxf=r"gff3|gtf",


rule reference_ensembl_annotation_download:
    """download genome gff"""
    output:
        temp(
            "<tmp>/reference_ensembl_annotation_download/{species}.{build}.{release}.{gxf}"
        ),
    log:
        "<log>/reference_ensembl_annotation_download/{species}.{build}.{release}.{gxf}.log",
    benchmark:
        "<benchmark>/curl/reference_ensembl_annotation_download_{species}.{build}.{release}.{gxf}.tsv"
    threads: 1
    params:
        species="{species}",
        build="{build}",
        release="{release}",
    wrapper:
        "v9.4.2/bio/reference/ensembl-annotation"


rule agat_default_config:
    """define agat configuration for gff"""
    output:
        temp("<tmp>/default_agat_config.yaml"),
    log:
        "<log>/agat_default_config.log",
    benchmark:
        "<benchmark>/agat/default_config.tsv"
    threads: 1
    shadow:
        "minimal"
    params:
        command="config",
        extra="",
    wrapper:
        "v9.6.0/bio/agat"


rule go_yq_update_gff_config:
    """define agat configuration for gff"""
    input:
        "<tmp>/default_agat_config.yaml",
    output:
        "<reference>/{species}.{build}.{release}/annotations/{species}.{build}.{release}.agat_gff3_config.yaml",
    log:
        "<log>/go_yq_update_config/{species}.{build}.{release}.gff3.log",
    benchmark:
        "<benchmark>/go_yq/update_config_{species}.{build}.{release}.gff3.tsv"
    threads: 1
    params:
        extra="--verbose",
        command="eval",
        expression=".cpu=6 | .log=false | .progress_bar=false | .minimum_chunk_size=1000",
    wrapper:
        "v9.14.0/utils/go-yq"


use rule go_yq_update_gff_config as go_yq_update_gtf_config with:
    """define agat configuration for gtf"""
    input:
        "<reference>/{species}.{build}.{release}/annotations/{species}.{build}.{release}.agat_gff3_config.yaml",
    output:
        "<reference>/{species}.{build}.{release}/annotations/{species}.{build}.{release}.agat_gtf_config.yaml",
    log:
        "<log>/go_yq_update_config/{species}.{build}.{release}.gtf.log",
    benchmark:
        "<benchmark>/go_yq/update_config_{species}.{build}.{release}.gtf.tsv"
    threads: 1
    params:
        extra="--verbose",
        command="eval",
        expression='.output_format = "GTF"',


use rule agat_default_config as agat_convert_sp_gxf2gxf with:
    """fix common ensembl format issues"""
    input:
        gff="<tmp>/reference_ensembl_annotation_download/{species}.{build}.{release}.{gxf}",
        config="<reference>/{species}.{build}.{release}/annotations/{species}.{build}.{release}.agat_{gxf}_config.yaml",
    output:
        out=temp("<tmp>/agat_convert_sp_gxf2gxf/{species}.{build}.{release}.{gxf}"),
    log:
        "<log>/agat_convert_sp_gxf2gxf/{species}.{build}.{release}.{gxf}.log",
    benchmark:
        "<benchmark>/agat/convert_sp_gxf2gxf_{species}.{build}.{release}.{gxf}.tsv"
    threads: 6
    shadow:
        "minimal"
    params:
        command="agat_convert_sp_gxf2gxf.pl",
        extra="",


use rule agat_default_config as agat_sq_filter_feature_from_fasta with:
    """ensure fasta and gff have the same contigs"""
    input:
        gff="<tmp>/agat_convert_sp_gxf2gxf/{species}.{build}.{release}.{gxf}",
        config="<reference>/{species}.{build}.{release}/annotations/{species}.{build}.{release}.agat_{gxf}_config.yaml",
        fasta="<genome_sequence>",
    output:
        o=temp(
            "<tmp>/agat_sq_filter_feature_from_fasta/{species}.{build}.{release}.{gxf}"
        ),
    log:
        "<log>/agat_sq_filter_feature_from_fasta/{species}.{build}.{release}.{gxf}.log",
    benchmark:
        "<benchmark>/agat/sq_filter_feature_from_fasta_{species}.{build}.{release}.{gxf}.tsv"
    threads: 6
    shadow:
        "minimal"
    params:
        extra="",
        command="agat_sq_filter_feature_from_fasta.pl",


use rule agat_default_config as agat_sp_filter_feature_by_attribute_value with:
    """filter out tsl na from genome annotation"""
    input:
        gff="<tmp>/agat_sq_filter_feature_from_fasta/{species}.{build}.{release}.{gxf}",
        config="<reference>/{species}.{build}.{release}/annotations/{species}.{build}.{release}.agat_{gxf}_config.yaml",
    output:
        o="<reference>/{species}.{build}.{release}/annotations/{species}.{build}.{release}.{gxf}",
    log:
        "<log>/agat_sp_filter_feature_by_attribute_value/{species}.{build}.{release}.{gxf}.log",
    benchmark:
        "<benchmark>/agat/sp_filter_feature_by_attribute_value_{species}.{build}.{release}.{gxf}.tsv"
    threads: 6
    shadow:
        "minimal"
    params:
        extra="--attribute 'transcript_support_level' --value 'NA' --test '='",
        command="agat_sp_filter_feature_by_attribute_value.pl",
