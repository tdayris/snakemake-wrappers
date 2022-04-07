#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Adrien Legendre
# SeqOIA
# 13/03/2020
# Parse statistic files from samtools idxstats, samtools flagstats, fastqc and snpeff and create a csv file to summarize all the statistics 
# Test with python3

import argparse
import os
import csv
import json
import re
import glob
import numpy as np
    
###############################################################################
#_______________FASTQC.TXT PARSING__________________

def fastqc_parsing(fastqc_file):
    fastqc_result = []
    with open(fastqc_file[0], 'r') as txt_file:
        fastqc = []
        fastqc = txt_file.readlines()
        fastqc[-1] = fastqc[-1].split('\n')[0]
        fastqc=fastqc[0].split(' ')
        for fastqc_f in fastqc:
            basic_stat=[]
            with open(fastqc_f,"r") as f:
                file_list=f.readlines()
                file_sublist=[]
                temp=[]
                for item in file_list: #__Sublist creation by spliting the file using the string >>END_MODULE
                    temp.append(item)
                    if item==">>END_MODULE\n":
                        file_sublist.append(temp)
                        temp=[]
                if temp:
                    file_sublist.append(temp)
                #Determine Sublist by parsing the text file with END_MODULE#####
                basic_stat=file_sublist[0]
                sample_result=[basic_stat[3].split('\t')[1][:-1], basic_stat[6].split('\t')[1][:-1], basic_stat[8].split('\t')[1][:-1],basic_stat[9].split('\t')[1][:-1]]
            fastqc_result.append(sample_result)
    return fastqc_result

#########################################################################################
#______________Flagstats file parsing_____________
def flagstats_parsing(list_i):
    flagstat_result = []
    with open(list_i[0], 'r') as txt_file:
        files = []
        files = txt_file.readlines()
        files[-1] = files[-1].split('\n')[0]
        files = files[0].split(' ')
    data_sample = []
    for v in files:
        data_file = []
        print("file :", v)
        with open(v,"r") as f:
            file_split = f.readlines()
            qc_pass = file_split[5].split(" ")[0]
            secondary = file_split[6].split(" ")[0]
            supplementary = file_split[7].split(" ")[0]
            duplicates = file_split[8].split(" ")[0]
            mapped = file_split[9].split(" ")[0]
            paired = file_split[10].split(" ")[0]
            read1 = file_split[11].split(" ")[0]
            read2 = file_split[12].split(" ")[0]
            properly = file_split[13].split(" ")[0]
            itself = file_split[14].split(" ")[0]
            singletons = file_split[15].split(" ")[0]
            mate_chr = file_split[16].split(" ")[0]
            mate_chr_5 = file_split[17].split(" ")[0]
            pct_mapped = float(mapped)*100/float(qc_pass)
            print("PCT :", pct_mapped)
            pct_properly = float(properly)*100/float(paired)
            pct_singletons = float(singletons)*100/float(qc_pass)
            data_file = [os.path.split(v)[1],qc_pass, secondary, supplementary, duplicates, mapped,pct_mapped, paired, read1, read2, properly, pct_properly, itself, singletons,pct_singletons, mate_chr, mate_chr_5 ]
        flagstat_result.append(data_file)
    return flagstat_result
###############################################################
#____________ Parser for all other QC files
def file_parsing(list_i, header, **kwargs):
    n_line = kwargs.get('n_line', 1)
    callable_loci = kwargs.get('callable_loci', False)
    index_line = []
    with open(list_i[0], 'r') as txt_file:
        files = []
        files = txt_file.readlines()
        files[-1] = files[-1].split('\n')[0]
        files = files[0].split(' ')
    with open(files[0], "r") as f:
        list_txt = f.read().splitlines()
        for k in range(len(list_txt)):
            lines = list_txt[k]
            if all(x in lines for x in header[1:]):
                for k in range(1, n_line + 1):
                    index_line.append(list_txt.index(lines) + k)
                lines = lines.split("\t")
                print("lines :",lines)
                print("Index :",list(lines.index(x) for x in header[1:]))
                index_result = list(lines.index(x) for x in header[1:])
    data = []
    for fi in files:
        with open(fi, "r") as f:
            print("FI:", fi)
            filename = os.path.split(fi)[1]
            list_txt = f.read().splitlines()
            if callable_loci == True:
                temp = [filename]
                for h in header:
                    for line in list_txt:
                        if h in line and "PCT_" in h :
                            temp.append(line.split(h.split("PCT_")[1] + " ")[-1])
                            print("temp :", temp)
                        if h in line and not "PCT_" in h:
                            temp.append(line.split(h+" ")[-1])
                            print("temp :", temp)
                data.append(temp)
            else:
                print("index_line :", index_line)
                for i in index_line:
                    temp = [filename]
                    lines = list_txt[i].split("\t")
                    print(temp, lines)
                    print("Index result :", index_result)
                    for k in index_result:
                        temp.append(lines[k])
                    print("temp :", temp)
                    data.append(temp)
    print("DATA :", data)
    return data
#################################################################
#______________ Calculation of total statistics _______________
def total_calculation(list_i, indice, chromosome_ratio, **kwargs):
    insert_size = kwargs.get('insert_size', False)
    total = ["Total"]
    print("LIST_I :", list_i)
    print("Calcul len total = ", list_i[0])
    print("len_total = ", len(list_i[0]))
    for k in range(1,len(list_i[0])):
        total.append(0)
    #check if all indices are present
    len_total = len(total) - 1 #because of the first element which doesn't count for total calculation --> not present in the json indice
    len_indice = len(indice["weighted_average"]) + len(indice["sum"]) + len(indice["average"]) + len(indice["string"])
    if len_total != len_indice:
        raise ValueError("Some indices are missing in the Indice json provided, data length = {}, indice length = {}".format(len_total, len_indice))
    #changement of value "?", ".",""
    for k in list_i:
        for i in range(1,len(k)):
            if k[i] in ["?", "-", ""]:
                k[i] = 0
            #### Int and float data type changement
            elif k[i] in ["FR", "TANDEM", "RF", "FIRST_OF_PAIR", "PAIR", "SECOND_OF_PAIR", "Merge_CORE-SPIKE_GRCh38"]:
                continue
            else:       #type(k[i]) not in [int,float]:
                try:
                    k[i] = float(k[i])
                except ValueError:
                    try:
                        k[i] = int(k[i])
                    except ValueError:
                        raise ValueError("{} can't be converted into float or int".format(k[i]))
    print("LIST_I_A:", list_i)
    #Calculation
    for x in range(1, len(total)):
        t = []
        i = []
        if x in indice["weighted_average"] and x not in indice["average"] + indice["sum"]:
            print("weighted:", x)
            for lines in list_i:
                print("lines :",lines[0])
                p = re.compile('_chr_(\d+|X|Y|MT)')
                print("regex :", re.findall(p, lines[0]))
                chrom = re.findall(p, lines[0])[0]
                print("T :",t)
                t.append(lines[x])
                print("i :", i)
                print("ratio :", chromosome_ratio[chrom])
                i.append(chromosome_ratio[chrom])
            total[x] = np.average(t, weights = i)
        elif x in indice["average"] and x not in indice["sum"] + indice["weighted_average"]:
            print("average:", x)
            if insert_size == True and x == len(total)-1:
                FR_count = 0
                RF_count = 0
                TANDEM_count = 0
                for lines in list_i:
                    if lines[x] == "FR":
                        FR_count += 1
                    if lines[x] == "RF_count":
                        RF_count += 1
                    if lines[x] == "TANDEM":
                        TANDEM_count += 1
                total[x] = "FR:{}/RF:{}/TANDEM:{}".format(FR_count,RF_count, TANDEM_count) 
            else:
                for lines in list_i:
                    t.append(lines[x])
                total[x] = np.average(t)
            #### besoin de traiter le cas "FR/TANDEM/FR"
        elif x in indice["sum"] and x not in indice["average"] + indice["weighted_average"]:
            print("sum:",x)
            for lines in list_i:
                total[x] += lines[x]
        elif x in indice["string"]:
            continue
        else:
            raise ValueError("Indice {} is in more than one category (average/weighted_average/sum)".format(i))
    print(total)
    return total
############################################################
#_____________ Ratio Chromosome preparation _______________
def chr_ratio_calculation(chromosome, f):
#dict creation from chr list
    chromosome_list = {}
    chromosome_ratio = {}
    for k in chromosome:
        chromosome_list[k] = 0
        chromosome_ratio[k] = 0
    # Check if file is valid ether .bed or .dict
    if not f.endswith(".dict") and  not f.endswith(".bed"):
        raise ValueError("{} is not a .bed or .dict file".format(f))
    with open(f, "r") as fi:
        data = fi.readlines()
        if f.endswith(".dict"):
            for lines in data:
                lines =  lines.split("\t")
                if "SN" in lines[1] and lines[1].split('SN:')[1] in chromosome_list:
                    chromosome_list[lines[1].split('SN:')[1]] = int(lines[2].split('LN:')[1])
        if f.endswith(".bed"):
            for lines in data:
                lines =  lines.split("\t")
                if lines[0] in chromosome_list:
                    chromosome_list[lines[0]] += int(lines[2]) - int(lines[1])
        # Calculation of chromosome ratio
        total_l = 0
        for k in chromosome_list.keys():
            total_l += chromosome_list[k]
        print("total length of target :", total_l)
        for k in chromosome_ratio.keys():
            chromosome_ratio[k] = (chromosome_list[k]/total_l)*100
    print("Chromosome ratio :",chromosome_ratio)
    return chromosome_ratio
############################################################
#_____________ Header of files statistiques ________________
wgs_header = ["Filename", "GENOME_TERRITORY", "MEAN_COVERAGE", "SD_COVERAGE", "MEDIAN_COVERAGE", "PCT_EXC_MAPQ", "PCT_EXC_DUPE", "PCT_EXC_UNPAIRED", "PCT_EXC_BASEQ","PCT_EXC_OVERLAP", "PCT_EXC_CAPPED", "PCT_EXC_TOTAL", "PCT_1X", "PCT_10X", "PCT_15X", "PCT_20X", "PCT_30X", "PCT_60X", "PCT_100X", "HET_SNP_SENSITIVITY", "HET_SNP_Q"]
indice_wgs = {"weighted_average": [2,5,6,7,8,9,10,11,12,13,14,15,16,17,18],
                  "average": [3,4,19,20],
                  "sum": [1],
                  "string":[]}


insert_size_header=["Filename", "MEDIAN_INSERT_SIZE", "MIN_INSERT_SIZE", "MAX_INSERT_SIZE", "MEAN_INSERT_SIZE", "STANDARD_DEVIATION", "READ_PAIRS", "PAIR_ORIENTATION"]
indice_insert_size = {"weighted_average": [],
                  "average": [1,2,3,4,5,7],
                  "sum": [6],
                  "string":[]}


doc_header = ["Filename", "total", "granular_median", "mean", "%_bases_above_1","%_bases_above_20", "%_bases_above_30","%_bases_above_60","%_bases_above_150", "%_bases_above_200"]
indice_doc = {"weighted_average": [3,4,5,6,7,8,9],
                  "average": [2],
                  "sum": [1],
                  "string":[]}


alignement_header = ["Filename","CATEGORY", "TOTAL_READS", "PF_READS", "PCT_PF_READS","PF_NOISE_READS", "PF_READS_ALIGNED", "PCT_PF_READS_ALIGNED","PF_ALIGNED_BASES", "PF_HQ_ALIGNED_READS", "PF_HQ_ALIGNED_BASES", "PF_HQ_ALIGNED_Q20_BASES", "PF_HQ_MEDIAN_MISMATCHES","PF_MISMATCH_RATE", "PF_HQ_ERROR_RATE", "PF_INDEL_RATE", "MEAN_READ_LENGTH", "READS_ALIGNED_IN_PAIRS", "PCT_READS_ALIGNED_IN_PAIRS", "PF_READS_IMPROPER_PAIRS", "PCT_PF_READS_IMPROPER_PAIRS", "BAD_CYCLES", "PCT_ADAPTER","STRAND_BALANCE", "PCT_CHIMERAS"]
indice_alignement = {"weighted_average": [4,7,18,20,23],
                  "average": [12,13,14,15,16,21,22,24],
                  "sum": [2,3,5,6,8,9,10,11,17,19],
                  "string":[1]}


callable_loci_header = ["Filename","REF_N","CALLABLE", "NO_COVERAGE", "LOW_COVERAGE", "EXCESSIVE_COVERAGE", "POOR_MAPPING_QUALITY"]
indice_callable_loci = {"weighted_average": [],
                  "average": [],
                  "sum": [1,2,3,4,5,6],
                  "string":[]}


markdup_header=["Filename", "UNPAIRED_READS_EXAMINED", "READ_PAIRS_EXAMINED", "SECONDARY_OR_SUPPLEMENTARY_RDS","UNMAPPED_READS", "UNPAIRED_READ_DUPLICATES", "READ_PAIR_DUPLICATES","READ_PAIR_OPTICAL_DUPLICATES", "PERCENT_DUPLICATION", "ESTIMATED_LIBRARY_SIZE"]
indice_markdup = {"weighted_average": [8],
                  "average": [7],
                  "sum": [1,2,3,4,5,6,9],
                  "string":[]}


quality_yield_header=["Filename","TOTAL_READS", "PF_READS", "READ_LENGTH", "TOTAL_BASES","PF_BASES", "Q20_BASES", "PF_Q20_BASES", "Q30_BASES", "PF_Q30_BASES", "Q20_EQUIVALENT_YIELD", "PF_Q20_EQUIVALENT_YIELD"]
indice_quality_yield = {"weighted_average": [],
                        "average": [3],
                        "sum": [1,2,4,5,6,7,8,9,10,11],
                        "string": []}


flagstats_header = ["Filename", "QC-passed reads", "Secondary", "Supplementary alignment", "Duplicates", "Mapped","PCT mapped", "Paired in Sequencing", "Read1", "Read2", "Properly Paired", "PCT properly paired", "Properly paired with itself and mate mapped", "Singletons", "PCT singletons", "Reads with mate mapped to a different chr", "Reads with mate mapped to a different chr (mapQ>=5)"]
indice_flagstats = {"weighted_average": [],
                    "average": [6,11,14],
                    "sum": [1,2,3,4,5,7,8,9,10,12,13,15,16],
                    "string":[]}

fastqc_header=["Filename","Total Sequences", "Sequence length", "%GC"]
indice_fastqc = {"weighted_average": [],
                 "average": [2,3],
                 "sum": [1],
                 "string": []}
targeted_pcr_header=["Filename","CUSTOM_AMPLICON_SET","AMPLICON_TERRITORY","ON_AMPLICON_BASES","NEAR_AMPLICON_BASES","OFF_AMPLICON_BASES","PCT_AMPLIFIED_BASES","PCT_OFF_AMPLICON","ON_AMPLICON_VS_SELECTED","MEAN_AMPLICON_COVERAGE","FOLD_ENRICHMENT","PF_SELECTED_PAIRS","PF_SELECTED_UNIQUE_PAIRS","ON_TARGET_FROM_PAIR_BASES","TARGET_TERRITORY","GENOME_SIZE","TOTAL_READS","PF_READS","PF_BASES","PF_UNIQUE_READS","PF_UQ_READS_ALIGNED","PF_BASES_ALIGNED","PF_UQ_BASES_ALIGNED","ON_TARGET_BASES","PCT_PF_READS","PCT_PF_UQ_READS","PCT_PF_UQ_READS_ALIGNED","MEAN_TARGET_COVERAGE","MEDIAN_TARGET_COVERAGE","MAX_TARGET_COVERAGE","MIN_TARGET_COVERAGE","ZERO_CVG_TARGETS_PCT","PCT_EXC_DUPE","PCT_EXC_ADAPTER","PCT_EXC_MAPQ","PCT_EXC_BASEQ","PCT_EXC_OVERLAP","PCT_EXC_OFF_TARGET","FOLD_80_BASE_PENALTY","PCT_TARGET_BASES_1X","PCT_TARGET_BASES_2X","PCT_TARGET_BASES_10X","PCT_TARGET_BASES_20X","PCT_TARGET_BASES_30X","PCT_TARGET_BASES_40X","PCT_TARGET_BASES_50X","PCT_TARGET_BASES_100X","AT_DROPOUT","GC_DROPOUT","HET_SNP_SENSITIVITY","HET_SNP_Q","SAMPLE","LIBRARY","READ_GROUP"]
indice_targeted = {"weighted_average": [],
                   "average": [6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53],
                   "sum": [2,3,4,5],
                   "string": [1]}

def main(output, chromosome, qc_list, f, wts_flag):
    chromosome_ratio = chr_ratio_calculation(chromosome,f)
    if wts_flag == False and wes_flag == False:
        qc_order = ["FASTQC", "Picard CollectInsertSizeMetrics statistics", "Samtools Flagstats", "GATK Depth Of Coverage (Only Autosome)", "Picard Collect Alignement Summary Metrics", "GATK Callable Loci Metrics", "PICARD Collect WGS Metrics (Only Autosome)", "PICARD MarkeDuplicates Metrics", "PICARD Collect Quality Yield Metrics"]
    elif wts_flag == True and wes_flag == False:
        qc_order = ["FASTQC", "Picard CollectInsertSizeMetrics statistics", "Samtools Flagstats", "Picard Collect Alignement Summary Metrics", "PICARD MarkeDuplicates Metrics", "PICARD Collect Quality Yield Metrics"]
    elif wts_flag == False and wes_flag == True:
        qc_order = ["FASTQC", "Picard CollectInsertSizeMetrics statistics", "Samtools Flagstats", "GATK Depth Of Coverage (Only Autosome)", "Picard Collect Alignement Summary Metrics", "GATK Callable Loci Metrics", "PICARD MarkeDuplicates Metrics", "PICARD Collect Quality Yield Metrics", "PICARD Targeted PCR Metrics"]
    with open(output, "w+") as f:
        writer=csv.writer(f, quoting=csv.QUOTE_ALL)
        for qc in qc_order:
            print("Json key :", qc)
            ##
            indice = qc_list[qc]["indice"]
            files = qc_list[qc]["file"]
            header = qc_list[qc]["header"]
            ##
            if qc == "FASTQC":
                data = fastqc_parsing(files)
                data.append(total_calculation(data,indice,chromosome_ratio))
            elif qc == "Samtools Flagstats":
                print("Samtools Flagstats")
                data = flagstats_parsing(files)
                data.append(total_calculation(data,indice,chromosome_ratio))
            elif qc == "GATK Callable Loci Metrics":
                print("GATK Callable Loci Metrics")
                data = file_parsing(files, header, callable_loci=True)
                data.append(total_calculation(data,indice,chromosome_ratio))
            elif qc == "Picard CollectInsertSizeMetrics statistics":
                print("Picard CollectInsertSizeMetrics statistics")
                data = file_parsing(files, header)
                data.append(total_calculation(data,indice,chromosome_ratio, insert_size=True))
            elif qc == "Picard Collect Alignement Summary Metrics":
                print("Picard Collect Alignement Summary Metrics")
                data = file_parsing(files, header, n_line=3)
                for k in ["FIRST_OF_PAIR", "SECOND_OF_PAIR", "PAIR"]:
                    temp = []
                    for lines in data:
                        if k in lines:
                            temp.append(lines)
                    total_temp = total_calculation(temp,indice,chromosome_ratio)
                    total_temp[1] = k
                    data.append(total_temp)
            else:
                data = file_parsing(files, header)
                data.append(total_calculation(data,indice,chromosome_ratio))
            writer.writerow([qc])
            writer.writerow(header)
            for lines in data:
                print("l :", lines)
                writer.writerow(lines)
            writer.writerow([])
    print("FINAL_CSV_PATH :", output)
    return output

##############################################################
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Parse statistic files from samtools idxstats, samtools flagstats, fastqc and snpeff and create a csv file to summarize all the statistics')
    parser.add_argument('-doc', '--depth_of_coverage', type=str, nargs=1,action='append', help='Path of File containing list of all gatk3 depth of coverage summary file', required=False)
    parser.add_argument('-f', '--fastqc', type=str, nargs=1,action='append', help='Path of File containing list of all fastqc file', required=True)
    parser.add_argument('-fl', '--flagstats', type=str, nargs=1,action='append', help='Path of File containing list of all flagstat files, must be a dict python type', required=True)
    parser.add_argument('-is', '--insert_size', type=str, nargs=1,action='append', help='Path of File containing list of all picard CollectInsertSizeMetrics file', required=True)
    parser.add_argument('-o', '--output', type=str, nargs=1, help='path for csv output file', required=True)
    parser.add_argument('-am', '--alignement_metrics', type=str, nargs=1,action='append', help='Path of File containing list of all picard alignement metrics file', required=True)
    parser.add_argument('-cl', '--callable_loci', type=str, nargs=1,action='append', help='Path of File containing list of all GATK callable loci file', required=False)
    parser.add_argument('-wgs', '--wgs_metrics', type=str, nargs=1,action='append', help='Path of File containing list of all picard collect wgs metrics file', required=False)
    parser.add_argument('-m', '--markedup_metrics', type=str, nargs=1,action='append', help='Path of File containing list of all picard markedup metrics file', required=True)
    parser.add_argument('-qy', '--quality_yield', type=str, nargs=1,action='append', help='Path of File containing list of all picard quality yield metrics file', required=True)
    parser.add_argument('-b', '--bed', type=str, nargs=1,action='append', help='Path of dict or bed file for the target', required=True)
    parser.add_argument('-chr', '--chromosome', type=str, nargs=1,action='append', help='list of chromosome', required=False, default = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y", "MT"])
    parser.add_argument('-wts', '--wts_analysis', type=bool, help='options for WTS analysis', required=False, default = False)
    parser.add_argument('-wes', '--wes_analysis', type=bool, help='options for WES analysis', required=False, default = False)
    parser.add_argument('-tp', '--targeted_pcr', type=str, nargs=1, action='append', help='options for WES True', required=False)
    ####################
    args = parser.parse_args()
    markdup=args.markedup_metrics[0]
    alignement=args.alignement_metrics[0]
    fastqc=args.fastqc[0]
    flagstats=args.flagstats[0]
    insert_size=args.insert_size[0]
    quality_yield=args.quality_yield[0]
    output=args.output[0]
    bed = args.bed[0][0]
    chromosome = args.chromosome
    wts_flag = args.wts_analysis
    wes_flag = args.wes_analysis
    if wts_flag == False and wes_flag == False:
        wgs=args.wgs_metrics[0]
        doc=args.depth_of_coverage[0]
        callable_loci=args.callable_loci[0]
    if wes_flag == True:
        callable_loci=args.callable_loci[0]
        doc=args.depth_of_coverage[0]
        targeted = args.targeted_pcr[0]
    ###
    if wts_flag == False and wes_flag == False:
        qc_list = {"FASTQC" : {"file": fastqc, "indice": indice_fastqc, "header": fastqc_header},
               "Picard CollectInsertSizeMetrics statistics": {"file": insert_size, "indice": indice_insert_size, "header": insert_size_header},
               "Samtools Flagstats" : {"file": flagstats, "indice": indice_flagstats, "header": flagstats_header},
               "GATK Depth Of Coverage (Only Autosome)": {"file": doc, "indice": indice_doc, "header": doc_header},
               "Picard Collect Alignement Summary Metrics": {"file": alignement, "indice": indice_alignement, "header": alignement_header},
               "GATK Callable Loci Metrics" : {"file": callable_loci, "indice": indice_callable_loci, "header": callable_loci_header},
               "PICARD Collect WGS Metrics (Only Autosome)": {"file": wgs, "indice": indice_wgs, "header": wgs_header},
               "PICARD MarkeDuplicates Metrics": {"file": markdup, "indice": indice_markdup, "header": markdup_header},
               "PICARD Collect Quality Yield Metrics": {"file": quality_yield, "indice": indice_quality_yield, "header": quality_yield_header}}
    elif wts_flag == True and wes_flag == False:
        qc_list = {"FASTQC" : {"file": fastqc, "indice": indice_fastqc, "header": fastqc_header},
               "Picard CollectInsertSizeMetrics statistics": {"file": insert_size, "indice": indice_insert_size, "header": insert_size_header},
               "Samtools Flagstats" : {"file": flagstats, "indice": indice_flagstats, "header": flagstats_header},
               "Picard Collect Alignement Summary Metrics": {"file": alignement, "indice": indice_alignement, "header": alignement_header},
               "PICARD MarkeDuplicates Metrics": {"file": markdup, "indice": indice_markdup, "header": markdup_header},
               "PICARD Collect Quality Yield Metrics": {"file": quality_yield, "indice": indice_quality_yield, "header": quality_yield_header}}
    elif  wts_flag == False and wes_flag == True:
         qc_list = {"FASTQC" : {"file": fastqc, "indice": indice_fastqc, "header": fastqc_header},
               "Picard CollectInsertSizeMetrics statistics": {"file": insert_size, "indice": indice_insert_size, "header": insert_size_header},
               "Samtools Flagstats" : {"file": flagstats, "indice": indice_flagstats, "header": flagstats_header},
               "GATK Depth Of Coverage (Only Autosome)": {"file": doc, "indice": indice_doc, "header": doc_header},
               "Picard Collect Alignement Summary Metrics": {"file": alignement, "indice": indice_alignement, "header": alignement_header},
               "GATK Callable Loci Metrics" : {"file": callable_loci, "indice": indice_callable_loci, "header": callable_loci_header},
               "PICARD MarkeDuplicates Metrics": {"file": markdup, "indice": indice_markdup, "header": markdup_header},
               "PICARD Collect Quality Yield Metrics": {"file": quality_yield, "indice": indice_quality_yield, "header": quality_yield_header},
               "PICARD Targeted PCR Metrics": {"file": targeted, "indice": indice_targeted, "header": targeted_pcr_header}}
    ##
    main(output, chromosome, qc_list, bed, wts_flag)
