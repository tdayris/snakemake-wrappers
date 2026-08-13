rule test_mashmap_reflist:
    input:
        ref="reference.txt",
        query="read.fasta.gz",
    output:
        "mashmap.out",
    log:
        "logs/mashmap.log",
    params:
        extra="-s 500 --pi 99",
    wrapper:
        "master/bio/mashmap"
