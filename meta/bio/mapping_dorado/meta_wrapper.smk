rule download_basecalling_model_for_dorado:
    output:
        model=temp(directory("<temp>/download_basecalling_model_for_dorado/{model}")),
    log:
        "<logs>/download_basecalling_model_for_dorado/{model}",
    benchmark:
        "<benchmarks>/download_basecalling_model_for_dorado/{model}.tsv"
    container:
        "docker://quay.io/sangerpathogens/cuda_dorado@sha256:f2a15f1caa4c1ded33c0774f0771a3930139e340ef3c4ade3c36f2bd002f1b18"
    threads: 1
    params:
        model="{model}",
        model_dir=lambda wildcards, output: subpath(output[0], parent=True),
    shell:
        "dorado download "
        "--model '{params.model}' "
        "--models-directory '{params.model_dir}' "
        "> {log:q} 2>&1"


rule basecalling_with_dorado:
    """call bases with dorado"""
    input:
        pod5="<temp>/fast5_to_pod5_conversion/<sample>.pod5",
        model="<temp>/download_basecalling_model_for_dorado/rna004_130bps_fast@v5.2.0",
    output:
        ubam=temp("<temp>/basecalling_with_dorado/<sample>.ubam"),
    log:
        "<logs>/basecalling_with_dorado/<sample>.log",
    benchmark:
        "<benchmarks>/basecalling_with_dorado/<sample>.tsv"
    container:
        "docker://quay.io/sangerpathogens/cuda_dorado@sha256:f2a15f1caa4c1ded33c0774f0771a3930139e340ef3c4ade3c36f2bd002f1b18"
    threads: 20
    params:
        model_dir=lambda wildcards, input: subpath(input.model, parent=True),
        model=lambda wildcards, input: subpath(input.model, basename=True),
        extra="--recursive --trim ",
    shell:
        "dorado basecaller {params.extra} "
        "--models-directory {params.model_dir:q} "
        "{params.model:q} {input.pod5:q} "
        "> {output.ubam:q} 2> {log}"


rule index_genome_fasta_sequence_with_minimap2:
    """in order to align reads through dorado"""
    input:
        target="<genome_sequence>",
        target_idx="<genome_sequence_index>",
    output:
        "<resources>/<species>.<build>.<release>/minimap2/<species>.<build>.<release>.<seqtype>.mmi",
    log:
        "<logs>/index_genome_fasta_sequence_with_minimap2/<species>.<build>.<release>.<seqtype>.log",
    benchmark:
        "<benchmarks>/minimap2_index/index_genome_fasta_sequence_with_minimap2_<species>.<build>.<release>.<seqtype>.tsv",
    threads: 20
    params:
        extra="",
    wrapper:
        "v9.9.0/bio/minimap2/index"


rule align_reads_with_dorado:
    """using minimap2 in the background"""
    input:
        ubam="<temp>/basecalling_with_dorado/{sample}.ubam",
        index="<resources>/{species}.{build}.{release}/minimap2/{species}.{build}.{release}.{seqtype}.mmi",
    output:
        bam=temp(
            "<temp>/align_reads_with_dorado/{sample}.{seqtype}.{species}.{build}.{release}.bam"
        ),
        summary=temp(
            "<temp>/align_reads_with_dorado/{sample}.{seqtype}.{species}.{build}.{release}.summary"
        ),
    log:
        "<logs>/align_reads_with_dorado/{sample}.{seqtype}.{species}.{build}.{release}.log",
    benchmark:
        "<benchmarks>/align_reads_with_dorado/{sample}.{seqtype}.{species}.{build}.{release}.tsv"
    container:
        "quay.io/sangerpathogens/cuda_dorado@sha256:f2a15f1caa4c1ded33c0774f0771a3930139e340ef3c4ade3c36f2bd002f1b18"
    threads: 20
    params:
        extra="--emit-summary ",
        outdir=lambda wildcards, output: str(
            f"{subpath(output.bam, parent= True)}/"
            f"{wildcards.sample}.{wildcards.seqtype}."
            f"{wildcards.species}.{wildcards.build}.{wildcards.release}"
        ),
    shell:
        "dorado aligner {params.extra} "
        "--threads {threads} "
        "--output-dir {params.outdir:q} "
        "{input.index:q} {input.ubam:q} > {log} 2>&1 && "
        'bam=$(find {params.outdir:q} -type f -name "*.bam") && '
        'summary=$(find {params.outdir:q} -type f -name "alignment_summary.txt") && '
        "mv --verbose ${{bam}} {output.bam:q} >> {log:q} 2>&1 && "
        "mv --verbose ${{summary}} {output.summary:q} >> {log:q} 2>&1 "


rule sort_aligned_reads_with_sambamba:
    """gain space and processing speed"""
    input:
        "<temp>/align_reads_with_dorado/<sample>.<seqtype>.<species>.<build>.<release>.bam",
    output:
        "<results>/<species>.<build>.<release>/Dorado.<seqtype>/<sample>.bam",
    log:
        "<logs>/sort_aligned_reads_with_sambamba/<sample>.<seqtype>.<species>.<build>.<release>.log",
    benchmark:
        "<benchmarks>/sambamba_sort/sort_aligned_reads_from_dorado_<sample>.<seqtype>.<species>.<build>.<release>.tsv",
    threads: 15
    params:
        extra="",
    wrapper:
        "v3.11.0/bio/sambamba/sort"


rule index_sorted_reads_with_sambamba:
    """gain processing speed"""
    input:
        "<results>/<species>.<build>.<release>/Dorado.<seqtype>/<sample>.bam",
    output:
        "<results>/<species>.<build>.<release>/Dorado.<seqtype>/<sample>.bam.bai",
    log:
        "<logs>/index_sorted_reads_with_sambamba/<sample>.<seqtype>.<species>.<build>.<release>.log",
    benchmark:
        "<benchmarks>/sambamba_index/index_sorted_reads_from_dorado.<sample>.<seqtype>.<species>.<build>.<release>.tsv",
    threads: 3
    params:
        extra="",
    wrapper:
        "v6.1.0/bio/sambamba/index"
        
