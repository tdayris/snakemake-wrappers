"""
grep_pass (optional) : grep -P "^#" {patient}_CTC_all.cutadapt.sorted.rmmarkdup.hc.snps.filtered.vcf >{patient}_CTC_all.cutadapt.sorted.rmmarkdup.hc.snps.filtered.PASS.vcf ; grep -P "\tPASS\t" {patient}_CTC_all.cutadapt.sorted.rmmarkdup.hc.snps.filtered.vcf >>{patient}_CTC_all.cutadapt.sorted.rmmarkdup.hc.snps.filtered.PASS.vcf
"""

rule grep_pass_mutect2_ctc:
    input:
        "gatk/mutect2/{sample}_{version}_{manip}_{nb}.vcf"
    output:
        temp("gatk/mutect2/{sample}_{version}_{manip}_{nb}.filtered.vcf"),
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        tmpdir=tmp,
        time_min=get_20min_per_attempt,
    params:
        '-P "^#|\tPASS\t"'
    log:    
        "logs/grep_pass_mutect2/{sample}_{version}_{manip}_{nb}.log"
    shell:
        "grep {params} {input} > {output} 2> {log}"


rule grep_pass_mutect2_wbc:
    input:
        "gatk/mutect2/{sample}_{version}_{manip}.vcf"
    output:
        temp("gatk/mutect2/{sample}_{version}_{manip}.filtered.vcf"),
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        tmpdir=tmp,
        time_min=get_20min_per_attempt,
    params:
        '-P "^#|\tPASS\t"'
    log:    
        "logs/grep_pass_mutect2/{sample}_{version}_{manip}.log"
    shell:
        "grep {params} {input} > {output} 2> {log}"