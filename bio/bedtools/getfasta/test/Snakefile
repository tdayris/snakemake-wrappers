rule test_bedtools_getfasta:
    input:
        ref="reference.fasta",
        bed="intervals.bed"  # This can be a VCF/BED/GFF
    output:
        "sequences.fasta"  # This can be a TSV/FASTA/BED
    params:
        extra=""  # Do not set `-tsv` or `-bedOut`
    log:
        "logs/outfasta.log"
    wrapper:
        "master/bio/bedtools/getfasta"