rule tximport:
    input:
        **unpack(get_tximport)
    output:
        txi="032.tximport/{genome_build}.{genome_release}/txi.{comparison}.RDS"
    