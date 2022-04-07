#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Adrien Legendre
# SeqOIA
# 19/06/19
# Parse statistic from QC.py and summarize it 

import os
import csv
import glob
import datetime
import argparse
from bs4 import BeautifulSoup

pwd=os.getcwd()

parser = argparse.ArgumentParser(description='Create the QC resume for a SeqOIA run with all the samples')
parser.add_argument('-i', '--sample_QC_files', type=str, nargs='*', action='append', help='QC files from QC.py for all samples', required=True)
parser.add_argument('-o', '--output', type=str, nargs=1,action='store', help='csv output file', required=True)
parser.add_argument('-s', '--snpeff', type=str, nargs='+',action='append', help='Path of File containing list of all snpeff statistic file', required=True)
parser.add_argument('-vm', '--variant_metrics', type=str, nargs='+',action='append', help='Path of File containing list of all picard collect variant metrics file', required=True)
parser.add_argument('-l', '--chromosome_length', type=str,action='append', default=[248956422, 242193529, 198295559, 190214555, 181538259, 170805979,159345973, 145138636, 138394717, 133797422, 135086622, 133275309, 114364328, 107043718, 101991189, 90338345, 83257441, 80373285, 58617616, 64444167, 46709983, 50818468, 156040895, 57227415, 16569], nargs='+', help='list of chromosome length to calculate depth', required=False)
parser.add_argument('-r', '--run', type=str, nargs = '+', action='append', help='Run ID', required=False)
parser.add_argument('-p', '--pipeline_version', type=str, nargs = 1, action='append', help='Pipeline version for QC versioning', required=True)
args = parser.parse_args()
snpeff=args.snpeff[0]
variant_metrics=args.variant_metrics[0]
sample_QC_files = args.sample_QC_files[0]
print("sample_QC_files :", sample_QC_files, type (sample_QC_files))
output = args.output[0]
chr_length = args.chromosome_length
version = args.pipeline_version[0]

if args.run:
    run = args.run[0]
    print("Run", run)
########## length of reference ################
total_chr_length=0
for k in range(len(chr_length)):
        total_chr_length = total_chr_length + int(chr_length[k])
############ ratio chromosome #################
chr_ratio = []
for k in chr_length:
        chr_ratio.append(k/total_chr_length)
########## ratio autosome ###############
length_autosome = 0
for k in range(23):
        length_autosome = length_autosome + int(chr_length[k])
autosome_ratio = []
for k in range(23):
        autosome_ratio.append(chr_length[k]/length_autosome)

#########################################################################################
#___________________Snpeff stat file parsing___________________

snpeff_header=["Samples","Number of variants processed","SNV", "Indels", "Ratio SNV/Indels"]
snpeff_result=[]
ts_tv_global = []
het_hom_global = []
for snpeff_file in snpeff:
	with open(snpeff_file,"r") as f:
		chromosome = int(snpeff_file.split('chr_')[-1].split('_')[0])
		variant_number=[]
		SNV=[]
		ins=[]
		deletion=[]
		#ts_tv_ratio=[]
		#het_hom=[]
		stat = []
		sample_list = []
		soup=BeautifulSoup(f, "html.parser")
		text=[strings for strings in soup.stripped_strings]
		for k in range(len(text)):
			if text[k].startswith('Number of variants processed'):
				variant_number.append(int(text[k+2].replace(',','')))
			if text[k].startswith('SNP'):
				SNV.append(int(text[k+1].replace(',','')))
			if text[k].startswith('INS'):
				ins.append(text[k+1])
			if text[k].startswith('DEL'):
				deletion.append(text[k+1])
			if text[k].startswith('All variants'):
				lines = text[k+1].split('\n')
				nb_sample = len(lines[0].split(","))-2
				for x in range(nb_sample):
					sample_list.append(lines[0].split(",")[x+1])
					stat.append(lines[-1].split(",")[x+1])
				ts_tv_ratio= [chromosome] + stat
		ts_tv_global.append(ts_tv_ratio)
		#Indels calculation and Ratio SNV/indels
		indels=[int(ins[0].replace(',',''))+int(deletion[0].replace(',',''))]
		if indels[0] == 0:
			ratio=[0]	
		else:	
			ratio=[float(SNV[0]/int(indels[0]))]
		snpeff_result.append([os.path.split(snpeff_file)[1]]+ variant_number + SNV + indels + ratio)
	###########################################################################################
	#________________Calculation of Total line for Snpeff statistics_____________
total_ts_tv = ["Total ratio Ti/Tv"]
print("nb_sample :", nb_sample)
for k in range(nb_sample):
	total_ts_tv.append(0)

nbr_chr = len(ts_tv_global)
print("ts_tv :", ts_tv_global)
for lines in ts_tv_global:
	for k in range(nb_sample):
		total_ts_tv[k+1] = total_ts_tv[k+1] + (float(lines[k+1])/nbr_chr)
for i in range(len(sample_list)):
	if "WGS" in sample_list[i] or "WES" in sample_list[i]:
		print(sample_list[i])
		sample_list[i] = sample_list[i].split('.variant')[0]
	else:
		sample_list[i] = str(sample_list[i].split('.variant')[0]) + "_WTS"

sample_list = ["-"] + sample_list
print(nbr_chr)

nbr_chr = len(snpeff_result)
total_snpeff=["Total/Mean Ratios",0,0,0,0]
for k in range(len(snpeff_result)):
	for i in range(1,len(snpeff_result[k])):
		total_snpeff[i]=total_snpeff[i]+snpeff_result[k][i]
total_snpeff[-1] = int(total_snpeff[2])/int(total_snpeff[3])
snpeff_result.append(total_snpeff)
#####################PICARD VARIANT METRICS #############################################
variant_sum_header = ["Analysis","Total SNPS","Num in DB SNP", "Novel SNPS", "Filtered SNPS", "PCT DBSNP","DBNSP TITV","Novel TITV", "Total indels","Novel indels","Filtered Indels","PCT DBSNP indels","Num in DB SNP indels","DBSNP ins del ratio", "Novel in s del ratio", "Total Multiallelic SNPS","Num in DB SNP multiallelic", "Total complex indels","Num in DB SNP complex indels", "SNP Reference bias", "Num Singletons"]
variant_detail_header = ["Samples","Het_hom_ratio"]

summary = []
detail = []
total_detail = []
for k in range(nb_sample):
	detail.append([])
	total_detail.append([])

for files in variant_metrics:
	with open(files, "r") as txt_f:
		analysis =files.split("/")[-1].split("_chr_")[0]
		chromosome = files.split("/")[-1].split("_chr_")[-1].split("_variant")[0]
		if files.split(".")[-1] == "variant_calling_detail_metrics":
			list_txt = txt_f.readlines()
			for k in range(len(list_txt)):
				lines = list_txt[k]
				lines = lines.split("\t")
				if lines[0].startswith("SAMPLE_ALIAS"):
					for x in range(nb_sample):
						lines=list_txt[k+x+1]
						lines=lines.split("\t")
						detail[x].append([lines[0], chromosome, lines[1]])
		else:
			data = []
			list_txt = txt_f.readlines()
			for k in range(len(list_txt)):
				lines = list_txt[k]
				lines = lines.split("\t")
				if lines[0].startswith("TOTAL_SNPS"):
					lines=list_txt[k+1]
					lines=lines.split("\t")
					for k in range(len(lines)):
						if k == len(lines)-1:
							data.append(lines[k][:-1])
						else:
							data.append(lines[k])
					data = [chromosome] + data
					summary.append(data)
for k in range(len(detail)):
	sample = detail[k][0][0]
	if "WES-T" in sample or "WGS-C" in sample:
		sample = detail[k][0][0].split('.variant')[0]
	else:
		sample = str(detail[k][0][0].split('.variant')[0]) + "_WTS"
	het_hom = 0
	for x in range(len(detail[k])):
		chromosome = detail[k][x][1]
		if chromosome == "MT":
			if detail[k][x][2] == "?":
				het_hom = float(het_hom) + 0*chr_ratio[-1]
			else:
				het_hom = float(het_hom) + float(detail[k][x][2])*chr_ratio[-1]
		elif chromosome == "Y":
			if detail[k][x][2] == "?":
				het_hom = float(het_hom) + 0*chr_ratio[-2]
			else:
				het_hom = float(het_hom) + float(detail[k][x][2])*chr_ratio[-2]
		elif chromosome == "X":
			if detail[k][x][2] == "?":
				het_hom = float(het_hom) + 0*chr_ratio[-3]
			else:
				het_hom = float(het_hom) + float(detail[k][x][2])*chr_ratio[-3]
		else:
			if detail[k][x][2] == "?":
				het_hom = float(het_hom) + 0*chr_ratio[int(chromosome)-1]
			else:
				het_hom = float(het_hom) + float(detail[k][x][2])*chr_ratio[int(chromosome)-1]
		print("SAMPLE :", sample)
		print("HET :", het_hom)
	total_detail[k] = [sample, het_hom]

	###### MEAN and total % CALCULATION #############################################
total_sum = [analysis,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
for lines in summary:
	for k in range(1,len(lines)):
		if lines[k] == "?":
			lines[k] = 0
		if k in [1,2,3,4,8,9,10,12,16,17,18,20]:
			total_sum[k] = int(total_sum[k]) + int(lines[k])
		else:
			total_sum[k] = float(total_sum[k]) + (float(lines[k])/len(summary))
total_sum[5] = float(total_sum[5]) * 100
total_sum[11] = float(total_sum[11]) * 100
########################################### QC PARSING #######################################################


data_fastqc = []
data_idxstat = []
data_insert_size = []
data_flagstat = []
data_gatk_doc = []
data_alignement = []
data_callable = []
data_markedup=[]
data_yield = []
data_wgs = []
data_targeted = []
for k in sample_QC_files:
	data = []
	filename = k.split('/')[-1]
	analysis = k.split('/')[-3]
	with open(k, "r") as qc_file:
		qc_file = csv.reader(qc_file)
		for row in qc_file:
			if row and row[0].startswith('Total'):
				new_filename = filename.split('_S')
				if new_filename[1] == "M":
					row[0] = new_filename[0] + "_" + new_filename[1]
				else:
					row[0] = new_filename[0]		
				data.append(row)
			if row and row[0] == "Weighted average":
				row[0] = new_filename[0]			
				data.append(row)
			if row and row[0] == "Average":
				row[0] = new_filename[0]
				data.append(row)
		if "WTS" in filename:
			data_fastqc.append(data[0])
			data_insert_size.append(data[1])
			data_flagstat.append(data[2])
			data_alignement.append(data[3])
			data_alignement.append(data[4])
			data_alignement.append(data[5])
			data_markedup.append(data[6])
			data_yield.append(data[7])
		elif "WES" in filename:
			data_fastqc.append(data[0])
			data_insert_size.append(data[1])
			data_flagstat.append(data[2])
			data_gatk_doc.append(data[3])
			data_alignement.append(data[4])
			data_alignement.append(data[5])
			data_alignement.append(data[6])
			data_callable.append(data[7])
			data_markedup.append(data[8])
			data_yield.append(data[9])
			data_targeted.append(data[10])
		else:
			data_fastqc.append(data[0])
			data_insert_size.append(data[1])
			data_flagstat.append(data[2])
			data_gatk_doc.append(data[3])
			data_alignement.append(data[4])
			data_alignement.append(data[5])
			data_alignement.append(data[6])
			data_callable.append(data[7])
			data_wgs.append(data[8])
			data_markedup.append(data[9])
			data_yield.append(data[10])
		#Sort of Picar Alignment summary metrics using 3 values (Pair, First of pair and Second of pair)
		sort_alignement = []
		for index in ["FIRST_OF_PAIR", "SECOND_OF_PAIR", "PAIR"]:
			for k in data_alignement:
				if k[1] == index:
					sort_alignement.append(k)
			sort_alignement.append([""])
#####################################################
############### FLAGSTAT ############################
for k in range(len(data_flagstat)):
	data_flagstat[k][6] = (float(data_flagstat[k][5])/float(data_flagstat[k][1])) * 100
	data_flagstat[k][11] = (float(data_flagstat[k][10])/float(data_flagstat[k][7])) * 100
	data_flagstat[k][14] = (float(data_flagstat[k][13])/float(data_flagstat[k][7])) * 100
############### DATA QUALITY YIELD
for k in range(len(data_yield)):
	data_yield[k] = data_yield[k][0:8] + [(float(data_yield[k][7])/float(data_yield[k][5]))*100] + data_yield[k][8:10] + [(float(data_yield[k][9])/float(data_yield[k][5]))*100] + data_yield[k][10:]
#####################################################
for k in range(len(sort_alignement)):
	if len(sort_alignement[k]) > 2:
		for i in [4, 7, 18, 20, 22, 24]:
			sort_alignement[k][i] = float(sort_alignement[k][i])*100
### DATA CALLABLE % ################################
new_data_callable = []
total_length_chr = []
## Calcul of chr length for the tool
for k in data_callable:
	t = 0
	for i in range(1,7):
		t = t + int(k[i].split('.')[0])
	total_length_chr.append(t)
## Calcul of %
percent = []
for k in range(len(data_callable)):
	data_p = []
	print("data :", data_callable[k])
	print("total :", total_length_chr[k])
	for i in range(1,len(data_callable[k])):
		print(i)
		data_p.append((int(data_callable[k][i].split(".")[0]) * 100)/ int(total_length_chr[k]))
	percent.append(data_p)

for k in range(len(data_callable)):
	print("OLD :", data_callable[k])
	print("percent :", percent)
	data_callable[k] = [data_callable[k][0], data_callable[k][1], percent[k][0], data_callable[k][2], percent[k][1], data_callable[k][3], percent[k][2], data_callable[k][4], percent[k][3], data_callable[k][5], percent[k][4], data_callable[k][6], percent[k][5]]
	print("NEW :", data_callable[k])
###########################################
############# WGS #########################
for k in range(len(data_wgs)):
	for i in range(5,19):
		data_wgs[k][i] = float(data_wgs[k][i]) * 100
###########################################
############ MARKDUP ######################
for k in range(len(data_markedup)):
	data_markedup[k][8] = float(data_markedup[k][8]) * 100
###########################################
########### TARGETED ######################
print("TARGETED PCR: ", data_targeted)
print("len targeted: ", len(data_targeted))
if len(data_targeted) != 0:
	for k in [6,7,17,18,19,20,21,22,24,25,26,32,33,34,35,36,37,39,40,41,42,43,44,45,46]:
		data_targeted[0][k] = float(data_targeted[0][k])*100
################################################
markedup_header = ["Filename", "Unpaired reads Examined", "Read Paired Examined", "Secondary or Supplementary reads","Unmapped reads", "Unpaired read duplicate", "Read pair duplicates","Read pair optical duplicates", "Percent_duplicates", "Estimated library size"]
fastqc_header = ["Filename","Total Sequences", "Sequence length", "%GC"]
insert_size_header = ["Filename", "Median Insert Size Value", "Min Insert Size Value", "Max Insert Size Value", "Mean Insert Size Value", "Standard Deviation", "Read Pairs", "Pair Orientation"]
flagstat_header = ["Filename", "QC-passed reads", "Secondary", "Supplementary alignment", "Duplicates", "Mapped","PCT mapped", "Paired in Sequencing", "Read1", "Read2", "Properly Paired", "PCT properly paired", "Properly paired with itself and mate mapped", "Singletons", "PCT singletons", "Reads with mate mapped to a different chr", "Reads with mate mapped to a different chr (mapQ>=5)"]
gatk_doc_header = ["Samples", "Total aligned bases", "Median", "Mean Coverage", "%_bases_above_1","%_bases_above_20", "%_bases_above_30","%_bases_above_60","%_bases_above_150", "%_bases_above_200"] 
illumina_header = ["Run", "Lane Number", "Total Bases", "PF Bases", "Total Reads", "PF Reads", "Total Clusters", "PF Clusters", "Mean Cluster per tile", "SD Clusters per tile", "Mean PCT PF Clusters per tile","SD PCT PF cluster per tile", "Mean PF Clusters per tile", "SD PF Clusters per tile"]
callableloci_header = ["Samples","REF_N","PCT_REF_N","Callable", "PCT_Callable","No Coverage", "PCT_No coverage", "Low coverage", "PCT_Low coverage", "Excessive coverage", "PCT_Excessive coverage", "Poor mapping quality", "PCT_Poor mapping quality"]
alignement_header = ["Samples","CATEGORY", "Total Read", "PF reads", "PCT_PF_reads","PF_noise_read", "PF_read_aligned", "PCT_PF_reads_aligned","PF_aligned_bases", "PF_HQ_aligned_read", "PF_HQ_aligned_bases", "PF_HQ_aligned_q20_bases", "PF_HQ_aligned_median","PF_mismatch_rate", "PF_HQ_error_rate", "PF_indel_rate", "Mean_read_length", "Read_aligned_in_pairs", "PCT_reads_aligned_in_pairs", "PF_read_improper", "PCT_PF_read_improper", "Bad_cycles", "PCT_adapter","Strand_balance", "PCT_chimeras"]
yield_header = ["Filename","Total read", "PF reads", "Read length", "Total bases","PF bases", "Q20 bases", "PF Q20 bases","PCT PF Q20 bases", "Q30 bases", "PF Q30 bases","PCT PF Q30 bases", "Q20 Equivalent yield", "PF Q20 equivalent yield"]
wgs_header = ["Samples", "Genome territory","Mean Coverage","SD coverage","Median coverage","PCT_EXC_MAPQ", "PCT_EXC_DUPE", "PCT_EXC_UNPAIRED","PCT_EXC_BASEQ", "PCT_EXC_OVERLAP","PCT_EXC_CAPPED", "PCT_EXC_TOTAL", "PCT_1X", "PCT_10X", "PCT_15X", "PCT_20X", "PCT_30X", "PCT_60X", "PCT_100X", "HET_SNP_SENSITIVITY", "HET_SNP_Q"]
targeted_header = ["Filename","CUSTOM_AMPLICON_SET","AMPLICON_TERRITORY","ON_AMPLICON_BASES","NEAR_AMPLICON_BASES","OFF_AMPLICON_BASES","PCT_AMPLIFIED_BASES","PCT_OFF_AMPLICON","ON_AMPLICON_VS_SELECTED","MEAN_AMPLICON_COVERAGE","FOLD_ENRICHMENT","PF_SELECTED_PAIRS","PF_SELECTED_UNIQUE_PAIRS","ON_TARGET_FROM_PAIR_BASES","TARGET_TERRITORY","GENOME_SIZE","TOTAL_READS","PF_READS","PF_BASES","PF_UNIQUE_READS","PF_UQ_READS_ALIGNED","PF_BASES_ALIGNED","PF_UQ_BASES_ALIGNED","ON_TARGET_BASES","PCT_PF_READS","PCT_PF_UQ_READS","PCT_PF_UQ_READS_ALIGNED","MEAN_TARGET_COVERAGE","MEDIAN_TARGET_COVERAGE","MAX_TARGET_COVERAGE","MIN_TARGET_COVERAGE","ZERO_CVG_TARGETS_PCT","PCT_EXC_DUPE","PCT_EXC_ADAPTER","PCT_EXC_MAPQ","PCT_EXC_BASEQ","PCT_EXC_OVERLAP","PCT_EXC_OFF_TARGET","FOLD_80_BASE_PENALTY","PCT_TARGET_BASES_1X","PCT_TARGET_BASES_2X","PCT_TARGET_BASES_10X","PCT_TARGET_BASES_20X","PCT_TARGET_BASES_30X","PCT_TARGET_BASES_40X","PCT_TARGET_BASES_50X","PCT_TARGET_BASES_100X","AT_DROPOUT","GC_DROPOUT","HET_SNP_SENSITIVITY","HET_SNP_Q","SAMPLE","LIBRARY","READ_GROUP"]

with open(output, "w+") as qc_global:
		writer = csv.writer(qc_global, quoting=csv.QUOTE_ALL)
		writer.writerow(["Run","Date", "Version"])
		writer.writerow([run, datetime.datetime.now(), version])
		writer.writerow("")
		writer.writerow(["FASTQC STATS"])
		writer.writerow(fastqc_header)
		for k in data_fastqc:
			writer.writerow(k)
		writer.writerow("")
		writer.writerow(["PICARD COLLECT INSERT SIZE METRICS"])
		writer.writerow(insert_size_header)
		for k in data_insert_size:
			writer.writerow(k)
		writer.writerow("")
		writer.writerow(["FLAGSTAT STAT"])
		writer.writerow(flagstat_header)
		for k in data_flagstat:
			writer.writerow(k)
		writer.writerow("")
		writer.writerow(["SNPEFF STAT", "Only Autosome"])
		writer.writerow(snpeff_header)
		for k in snpeff_result:
			writer.writerow(k)
		writer.writerow("")
		writer.writerow(sample_list)
		writer.writerow(total_ts_tv)
		writer.writerow("")
		writer.writerow(["GATK DEPTH OF COVERAGE", "Only Autosome/WES Analysis on bed"])
		writer.writerow(gatk_doc_header)
		for k in data_gatk_doc:
			writer.writerow(k)
		writer.writerow("")
		writer.writerow(["PICARD ALIGNMENT METRICS"])
		writer.writerow(alignement_header)
		for k in sort_alignement:
			writer.writerow(k)
		writer.writerow("")
		writer.writerow(["GATK CALLABLE LOCI", "Only Autosome"])
		writer.writerow(callableloci_header)
		for k in data_callable:
			writer.writerow(k)
		writer.writerow("")
		writer.writerow(["PICARD WGS METRICS"])
		writer.writerow(wgs_header)
		for k in data_wgs:
			writer.writerow(k)
		writer.writerow("")
		writer.writerow(["PICARD VARIANTS METRICS"])
		writer.writerow(variant_sum_header)
		writer.writerow(total_sum)
		writer.writerow("")
		writer.writerow(variant_detail_header)
		for k in total_detail:
			writer.writerow(k)
		writer.writerow("")
		writer.writerow(["PICARD MARKEDUPLICATES METRICS"])
		writer.writerow(markedup_header)
		for k in data_markedup:
			writer.writerow(k)
		writer.writerow("")
		writer.writerow(["PICARD Collect Quality Yield Metrics"])
		writer.writerow(yield_header)
		for k in data_yield:
			writer.writerow(k)
		writer.writerow("")
		if len(data_targeted) != 0:
			writer.writerow(["PICARD Targeted PCR Metrics"])
			writer.writerow(targeted_header)
			for k in data_targeted:
				writer.writerow(k)
