rule agat_extract_cdna_sequences:
    input:
        gff="<resources>/annotation/{species}.{build}.{release}.dna.gtf",
        fasta="<resources>/sequence/{species}.{build}.{release}.dna.fa",
    output:
        o="<resources>/sequence/{species}.{build}.{release}.cdna.fa",
    log:
        "<logs>/agat_extract_transcripts_sequences/{species}.{build}.{release}.log",
    params:
        command="agat_sp_extract_sequences.pl",
        extra="-t cds",
    wrapper:
        "v7.3.0/bio/agat"


rule agat_remove_contigs_without_sequence:
    input:
        gff="<resources>/annotation/{species}.{build}.{release}.filtered.gtf",
        fasta="<resources>/sequence/{species}.{build}.{release}.dna.fa",
    output:
        o="<resources>/annotation/{species}.{build}.{release}.dna.gtf",
    log:
        "<logs>/agat_remove_contigs_without_sequence/{species}.{build}.{release}.log",
    params:
        command="agat_sq_filter_feature_from_fasta.pl",
        extra="--verbose",
    wrapper:
        "v7.3.0/bio/agat"


rule agat_remove_transcript_without_tsl:
    input:
        gff="<resources>/annotation/{species}.{build}.{release}.formatted.gtf",
    output:
        o=temp("<resources>/annotation/{species}.{build}.{release}.filtered.gtf"),
    params:
        command="agat_sp_filter_feature_by_attribute_value.pl",
        extra="--attribute 'transcript_support_level' --value '\"NA\"' --test '='",
    log:
        "<logs>/remove_transcript_without_tsl/{species}.{build}.{release}.log",
    wrapper:
        "v7.3.0/bio/agat"


rule agat_clean_gtf_common_format_errors:
    input:
        gff="<resources>/annotation/{species}.{build}.{release}.raw.gtf",
    output:
        out=temp("<resources>/annotation/{species}.{build}.{release}.formatted.gtf"),
    params:
        command="agat_convert_sp_gxf2gxf.pl",
        extra="",
    log:
        "<logs>/clean_gtf_common_format_errors/{species}.{build}.{release}.log",
    wrapper:
        "v7.3.0/bio/agat"


rule get_reference_genome_annotation:
    output:
        temp("<resources>/annotation/{species}.{build}.{release}.raw.gtf"),
    params:
        species="{species}",
        release="{release}",
        build="{build}",
    log:
        "<logs>/get_genome_annotation/{species}.{build}.{release}.log",
    cache: "omit-software"
    wrapper:
        "v7.4.0/bio/reference/ensembl-annotation"


rule get_reference_genome_sequence:
    output:
        "<resources>/sequence/{species}.{build}.{release}.dna.fa",
    params:
        species="{species}",
        datatype="dna",
        build="{build}",
        release="{release}",
    log:
        "<logs>/get_genome_sequence/.{species}.{build}.{release}.dna.log",
    cache: "omit-software"
    wrapper:
        "v5.10.0/bio/reference/ensembl-sequence"
