#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Adrien Legendre 
# 2020/11/20
# python3
import json
import csv
import argparse
import errno
import fcntl
import time

def identito(identito_csv, pid_crc, date, id_anon, run, output, analyse, sample_sheet):
	with open(sample_sheet, "r") as f:
		while True:
			try:
				fcntl.flock(f, fcntl.LOCK_EX| fcntl.LOCK_NB)
				break
			except IOError as e:
				if e.errno != errno.EAGAIN:
					raise
				else:
					time.sleep(1)
		for lines in f:
			lines = lines.split(",")
			if lines[1] == id_anon:
				id_labo = lines[0]
	with open(identito_csv, "r") as f:
		while True:
			try:
				fcntl.flock(f, fcntl.LOCK_EX| fcntl.LOCK_NB)
				break
			except IOError as e:
				if e.errno != errno.EAGAIN:
					raise
				else:
					time.sleep(1)
		data = f.readlines()
		identito = {}
		for k in range(1,len(data)):
			lines = data[k].split('","')
			if len(lines) > 2 and "rs" in lines[2]:
				identito[lines[2]] = {"Chromosome" : lines[0].split('"')[1], "Position": lines[1], "Call Information": lines[3], "Confidence": lines[4], "Confidence_SNP": lines[5], "Statut Puce": lines[6], "Statut NGS": lines[7], "Frequence NGS": lines[8], "Profondeur": lines[9], "Concordance": lines[10].split('"\n')[0]}
			if id_anon in lines[0]:
				nb_c = lines[1]
				nb_i = lines[2]
				concordance = lines[3]
		identito["Informations generales"] = {}
		identito["Informations generales"]["Prescription"] = pid_crc
		identito["Informations generales"]["Analyse"] = analyse
		identito["Informations generales"]["Date Analyse"] = date
		identito["Informations generales"]["Sample"] = {}
		identito["Informations generales"]["Sample"][id_labo] = {}
		identito["Informations generales"]["Sample"][id_labo]["ID"] = id_anon
		identito["Informations generales"]["Sample"][id_labo]["Run"] = run
		identito["Resultats"] = {"Nombre de Concordant": nb_c, "Nombre inconcordant": nb_i, "% Concordance": concordance}
	with open(output, "w+") as o:
		while True:
			try:
				fcntl.flock(o, fcntl.LOCK_EX| fcntl.LOCK_NB)
				break
			except IOError as e:
				if e.errno != errno.EAGAIN:
					raise
				else:
					time.sleep(1)
		json.dump(identito,o, ensure_ascii=False, indent=4, sort_keys = True)
	return output
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Python script which convert an identito csv file into a json')
	parser.add_argument('-i', '--identito', type=str, nargs=1,help='Path of identito csv', required=True)
	parser.add_argument('-p', '--pid_dossier', type=str, nargs=1,help='pid_crc of the prescription', required=True)
	parser.add_argument('-d', '--date_analysis', type=str, nargs=1,help='date of bioinformatic analysis', required=True)
	parser.add_argument('-id', '--sample_id', type=str, nargs=1,help='id_anon of the sample', required=True)
	parser.add_argument('-r', '--run_id', type=str, nargs=1,help='id of the sample run', required=True)
	parser.add_argument('-a', '--analysis', type=str, nargs=1,help='analysis id', required=True)
	parser.add_argument('-o', '--output_path', type=str, nargs=1,help='output_file', required=True)
	parser.add_argument('-c', '--csv', type=str, nargs=1,help='sample sheet path', required=True)
	args = parser.parse_args()
	##
	identito_csv = args.identito[0]
	print("identito:", identito_csv)
	sample_sheet = args.csv[0]
	pid_crc = args.pid_dossier[0]
	date = args.date_analysis[0]
	id_anon = args.sample_id[0].split("_S")[0]
	run = args.run_id[0]
	output = args.output_path[0]
	analyse = args.analysis[0]
	##
	identito(identito_csv, pid_crc, date, id_anon, run, output, analyse, sample_sheet)
