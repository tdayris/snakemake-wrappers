#!/usr/bin/env pythoni3
# -*- coding: utf-8 -*-

# Adrien Legendre
# SeqOIA
# 21/06/19
# Allow to create the identito-vigilance file for SeqOIA runs

import glob
import csv
import os
import argparse

parser = argparse.ArgumentParser(description='Create the identito-vigilance resume for SeqOIA runs')
parser.add_argument('-i', '--sample_identito', type=str, nargs=1, action='store', help='Unique file with all identito-vigilance informations for all the sample created using grep "SAMPLE" in the identito-vigilance file of the run ', required=True)
parser.add_argument('-v', '--checkvariant_file', type=str, nargs='*',action='append', help='File sample.checkvariant for all chromosomes', required=True)
parser.add_argument('-o', '--output', type=str, nargs=1,action='store', help='output file', required=True)
parser.add_argument('-tg', '--threshold_confidence_genotype', type=float, nargs=1,action='store', help='Threshold value for genotype confidence', required=True)
parser.add_argument('-ts', '--threshold_confidence_SNP', type=float, nargs=1,action='store', help='Threshold value for SNP confidence', required=True)
parser.add_argument('-r', '--rsID_position', type=str, nargs=1,action='store', help='file with all rsID and position', required=True)
parser.add_argument('-t', '--threshold_validation_ngs', type=int, action='store', help='threshold to validate NGS status', required=False, default = 80)
args = parser.parse_args()

identito = args.sample_identito[0]
variant = args.checkvariant_file[0]
output = args.output[0]
rsID_pos = args.rsID_position[0]
seuilNGS = args.threshold_validation_ngs
if args.threshold_confidence_genotype:
	confidenceTest = args.threshold_confidence_genotype[0]
else:
	confidenceTest = args.threshold_confidence_genotype
confidenceTestSNP = args.threshold_confidence_SNP[0]

split_path = variant[0].split('/')[-1].split('_S')
print("split :", split_path)
if variant[0].split('/')[-1].split('_')[1] == "M":
	sample = ('_').join([split_path[0], split_path[1]])
else:
	sample = split_path[0]

#check of confidence for all rsID in the file per sample, if the confidence is under the threshold then Statut and call info = No call
#Allow us to build a list with all rsID to find the position in rsID_pos
sample_rsID = []
error_conf = []
with open(identito, "r") as txt_file:
	txt_file = txt_file.readlines()
	variant_to_find = len(txt_file)
	print("Variant_to_find:", variant_to_find)
	for lines in txt_file:
		lines = lines.split(';')
		rs = lines[1]
		ref = lines[2]
		alt = lines[3]
		callinfo = lines[9]
		statut = lines[8]
		confidence = lines[7].replace(',','.')
		confidenceSNP = lines[12].replace(',','.')
		if not float(confidence) >= confidenceTest and float(confidenceSNP) >= confidenceTestSNP:
			callinfo = "No Call"
			statut = "No Call"
			error_conf.append([rs, ref, alt, callinfo, statut, confidence, confidenceSNP])
		sample_rsID.append([rs, ref, alt, callinfo, statut, confidence, confidenceSNP])
		print([rs, ref, alt, callinfo, statut, confidence, confidenceSNP])
	print("Sample_rsID: len=", len(sample_rsID))
print("ERROR CONFIDENCE VALUE AND SNP VALUE")
print(len(error_conf))
for k in error_conf:
	print(k)



#Extraction of chromosome, ref, alt  and postion for a specifique rsID in sample_rsID in the rsID_pos file
with open(rsID_pos, "r") as txt_file:
	txt_file = txt_file.readlines()	
	test = []
	for lines in sample_rsID:
		for lines_rsID_pos in txt_file:
			#if not working try to erase the [:-1]
			if lines[0] == lines_rsID_pos.split('\t')[2][:-1]:
				chromosome = lines_rsID_pos.split('\t')[0]
				position = lines_rsID_pos.split('\t')[1]
				lines.append([chromosome, int(position)+1])
				#lines.append([chromosome, position])
			else:
				continue
#Calcule of depth and ratio for rsID
new_sample_rsID = []
for path in variant:
	with open(path, "r") as checkvariant:
		checkvariant = checkvariant.readlines()
		chromosome = path.split('_')[-1].split('.')[0] 
		for lines in sample_rsID:
			for variants in checkvariant:
				variants = variants.split('\t')
				if variants[2] == str(lines[-1][1]) and chromosome == lines[-1][0]:
					depth = int(variants[4])
					if depth == 0:
						refRatio = 0
						altRatio = 0
					else:
						ref = lines[1]
						alt = lines[2]
						if ref == "A":
							refRatio = (int(variants[5])/depth)*100
						if ref == "T":
							refRatio = (int(variants[6])/depth)*100
						if ref == "C":
							refRatio = (int(variants[7])/depth)*100
						if ref == "G":
							refRatio = (int(variants[8])/depth)*100
						if alt == "A":
							altRatio = (int(variants[5])/depth)*100
						if alt == "T":
							altRatio = (int(variants[6])/depth)*100
						if alt == "C":
							altRatio = (int(variants[7])/depth)*100
						if alt == "G":
							altRatio = (int(variants[8])/depth)*100
					if refRatio >= seuilNGS:
						new_line = lines +["XX", refRatio, depth]
					elif altRatio >= seuilNGS:
						new_line = lines + ["YY", altRatio, depth]
					else:
						new_line = lines + ["XY", altRatio, depth]
					new_sample_rsID.append(new_line)
					print(new_line)
				else:
					continue
print("new_sample_rsID: len", len(new_sample_rsID))
concordant = 0
inconcordant = 0
nNocall = 0

for lines in new_sample_rsID:
	print("lines:", lines)
	if lines[4]!= "" and lines[-3] != "" and lines[4]!= "No Call" and lines[-3] != "No Call" and lines[4]!= "Invalid" and lines[-3] != "Invalid":
		if lines[4] == lines[-3]:
			concordant = concordant + 1
			lines.append("YES")
			print("CONCORDANT")
		else:
			inconcordant = inconcordant + 1
			lines.append("NO")
			print("NON_CONCORDANT")
	else:
		lines.append("No Call")
NoCall_NGS = variant_to_find - len(new_sample_rsID)
NoCall_puce = len(new_sample_rsID) - (concordant + inconcordant)
NoCall = NoCall_NGS + NoCall_puce
if concordant == 0 and inconcordant == 0:
	homologie = 0
else:
	homologie = (concordant/(concordant + inconcordant))*100
data_resume = [sample, concordant, inconcordant, NoCall, homologie]
print("data_resume:", data_resume)

with open(output, "w+") as output_file:
	writer=csv.writer(output_file, quoting=csv.QUOTE_ALL)
	writer.writerow(["Chromosome","Position","rsID","Call Information","Confidence","Confidence SNP","Statut Puce","Statut NGS","Frequence NGS","Profondeur","Concordance"])
	for k in new_sample_rsID:
		writer.writerow([k[7][0],k[7][1],k[0],k[4],k[5],k[6],k[3],k[8],k[9],k[10],k[-1]])
	writer.writerow("")
	writer.writerow(["Nombre de variants non retrouves dans le vcf:",NoCall_NGS])
	writer.writerow("")
	writer.writerow(data_resume)
