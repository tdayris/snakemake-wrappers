#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Adrien Legendre 
# 2020/08/27
#python3
import json
import csv
import argparse

def json_creation(input_f, output_f):
	with open(input_f, "r") as f:
		data = f.readlines()
		identito = {}
		for k in range(1,len(data)):
			lines = data[k].split('","')
			if len(lines) > 2 and "rs" in lines[2]:
				identito[lines[2]]= {"Chromosome" : lines[0].split('"')[1], "Position": lines[1], "Call Information": lines[3], "Confidence": lines[4], "Confidence_SNP": lines[5], "Statut Puce": lines[6], "Statut NGS": lines[7], "Frequence NGS": lines[8], "Profondeur": lines[9], "Concordance": lines[10].split('"\n')[0]}
		identito["Informations generales"] = {}
		identito["Informations generales"]["Prescription"] = "BCA4"
		identito["Informations generales"]["Analyse"] = "A00651_0034_WGS_trio_NA19238_19062020"
		identito["Informations generales"]["Date Analyse"] = "19062020"
		identito["Informations generales"]["Sample"] = "NA19238"
		identito["Resultats"] = {"Nombre de Concordant": "93", "Nombre inconcordant": "0", "% Concordance": "100%"}
	with open("/home/alegendre/710_WES-T_S7_identito.json", "w+") as f:
		json.dump(identito,f, ensure_ascii=False, indent=4, sort_keys = True)
	return output_f
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='VDM transcript analyse')
	parser.add_argument('-csv', '--sample_sheet', type=str, nargs=1,action='append', help='Path of sample sheet', required=True)
	parser.add_argument('-o', '--output', type=str, nargs=1,action='append', help='Path of output repertory', required=True)
	args = parser.parse_args()
	input_f = args.sample_sheet[0][0]
	output_f = args.output[0][0]
	#######
	json_creation(input_f, output_f)

