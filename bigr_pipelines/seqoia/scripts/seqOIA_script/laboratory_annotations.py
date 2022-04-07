# -*- coding: utf-8 -*-
# Adrien Legendre 
# 2019/08/26
from __future__ import unicode_literals
import os
import argparse
import sys
import csv
import json
import errno
import fcntl

def laboratory_name_extraction(sample_sheet):
	data = []
	duo_name = []
	with open(sample_sheet, "r") as csv_file:
		csv_data= csv.reader(csv_file)
		for row in csv_data:
			data.append(row)
	data = data[2:]
	for k in data:
		duo_name.append([k[0], k[1]])
	return duo_name

def csv_modification(qc_final, duo_name, sample_run):
	data = []
	with open(qc_final, "r") as csv_file:
		csv_list = csv.reader(csv_file)
		for row in csv_list:
			data.append(row)
	validation = [["-", "WGS-C", "-", "-", "-"], ["-", "WGS-T", "-", "-", "-"], ["-", "WES-T", "-", "-", "-"], ["-", "WTS", "-", "-", "-"]]
	for k in range(len(data)):
		nb_lines = len(data)
		if len(data[k]) > 0:
			if data[k][0] == "GATK DEPTH OF COVERAGE":
				i = 2
				while k+i < nb_lines and len(data[k+i]) != 0:
					if "WES-T" in data[k+i][0]:
						validation[2][2] = float(data[k+i][3])
						validation[2][1] = data[k+i][0]
					elif "WGS-T" in data[k+i][0]:
						validation[1][2] = float(data[k+i][3])
						validation[1][1] = data[k+i][0]
					elif "WGS-C" in data[k+i][0]:
						validation[0][2] = float(data[k+i][3])
						validation[0][1] = data[k+i][0]
					i = i + 1
			if data[k][0] == "PICARD WGS METRICS":
				i = 2
				while k+i < nb_lines and len(data[k+i]) != 0:
					if "WGS-C" in data[k+i][0]:
						validation[0][4] = float(data[k+i][14])
					elif "WGS-T" in data[k+i][0]:
						validation[1][4] = float(data[k+i][16])
					i = i + 1
			if data[k][0] == "PICARD Collect Quality Yield Metrics":
				i = 2
				while k+i < nb_lines and len(data[k+i]) != 0:
					if "WGS-T" in data[k+i][0]:
						validation[1][3] = float(data[k+i][10])
					elif "WGS-C" in data[k+i][0]:
						validation[0][3] = float(data[k+i][10])
					elif "WTS" in data[k+i][0]:
						validation[3][3] = float(data[k+i][10])
						validation[3][1] = data[k+i][0]
					i = i + 1
	for k in validation:
		print("validation:", k)
	for k in data:
		if len(k) != 0:
			if k[0] and k[0] in ["Samples","Filename"] and k[1] != "Number of variants processed":
				index = data.index(k)
				k = ["Sample Name"] + k
				data[index] = k
	for x in duo_name:
		sample_id = x[0]
		sample_name = x[1]
		for y in data:
			if len(y) != 0:
				if sample_name == y[0]:
					index = data.index(y)
					y = [sample_id] + y
					data[index] = y
	samples = []
	for x in duo_name:
		for k in range(len(validation)):
			if validation[k][1] == x[1]:
				validation[k][0] = x[0]
				samples.append(x)
###############################
	with open(output, "w") as csv_output:
		valid = []
		qc_valid = {}
		writer = csv.writer(csv_output, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
		for k in range(len(data)):
			data[1][0] = [sample_run]
			if k == 1:
				writer.writerow(data[k])
				writer.writerow([""])
				sample_done = []
				writer.writerow(["Sample Name", "Samples", "Profondeur moyenne (threshold: 30/60/150X)","Nombres de bases >= Q30 (threshold: 85/170/12.8Gb)", "Couverture mapq >=20 (threshold: 85%-15X/30X)"])
				for k in validation:
					writer.writerow(k)
					qc_valid[str(k[1])] = {"MEAN_COV": str(k[2]), "Q30": str(k[3]), "COUVERTURE" : str(k[4])}
			else:
				print("data :", data[k])
				writer.writerow(data[k])
	print("QC_VALID :", qc_valid)
	return output, qc_valid

def config_file_modification(config, qc_valid):
	with open(config, "r", encoding = 'utf-8') as f:
		while True:
			try:
				fcntl.flock(f, fcntl.LOCK_EX| fcntl.LOCK_NB)
				break
			except IOError as e:
				if e.errno != errno.EAGAIN:
					raise
				else:
					time.sleep(1)
		data = json.load(f)
		data["Compte_rendu"]["QC"] = qc_valid
		print(data["Compte_rendu"]["QC"])
	with open(config, "w", encoding = 'utf-8') as f:
		while True:
			try:
				fcntl.flock(f, fcntl.LOCK_EX| fcntl.LOCK_NB)
				break
			except IOError as e:
				if e.errno != errno.EAGAIN:
					raise
				else:
					time.sleep(1)
		json.dump(data,f, ensure_ascii=False, indent=4, sort_keys = True)
	return config

def sample_run_assoc(config):
	with open(config, "r", encoding = 'utf-8') as f:
		while True:
			try:
				fcntl.flock(f, fcntl.LOCK_EX| fcntl.LOCK_NB)
				break
			except IOError as e:
				if e.errno != errno.EAGAIN:
					raise
				else:
					time.sleep(1)
		data = json.load(f)
		sample_run = []
		for k in data["general_informations"]["SAMPLES"]:
			s = k.split("_")[1]
			sample_run.append(str(data["general_informations"][s]["RUN"]) + " : " + str(k))
	return sample_run
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='This python script add the Sample name to the QC final csv')
	parser.add_argument('-q', '--qc_final', type=str, nargs=1,help='Path of QC final file', required=True)
	parser.add_argument('-s', '--sample_sheet', type=str, nargs='+',help='Path of the sample sheet run', required=True)
	parser.add_argument('-o', '--output',type=str, nargs=1,help='Path of the output file', required=True)
	parser.add_argument('-c', '--config_file', type=str, nargs=1,help='Path of the output file', required=True)

	args = parser.parse_args()
	qc_final = args.qc_final[0]
	sample_sheet = args.sample_sheet
	output = args.output[0]
	config = args.config_file[0]
	duo_name = []
	for k in sample_sheet:
		duo_name = duo_name + laboratory_name_extraction(k)
	print("duo_name", duo_name)
	sample_run = sample_run_assoc(config)
	csv_modification(qc_final, duo_name, sample_run)
	qc_valid = csv_modification(qc_final, duo_name, sample_run)[1]
	config_file_modification(config, qc_valid)
