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

def qc_json(qc, pid_crc, date, json_pipe, output, analyse, pipe_version):
	with open(json_pipe, "r") as f:
		data = json.load(f)
		samples = {}
		samples["ID"] = {}
		samples["RUN"] = {}
		for k in data["general_informations"]["SAMPLES"]:
			print(k)
			samples["ID"][k.split('_S')[0]] = {}
			samples["ID"][k.split('_S')[0]]["RUN"] = data["general_informations"][k.split("_")[1]]["RUN"]
			samples["RUN"][data["general_informations"][k.split("_")[1]]["RUN"]] = {}
			samples["RUN"][data["general_informations"][k.split("_")[1]]["RUN"]]["LANE"] = data["general_informations"][k.split("_")[1]]["FLOWCELL"]["LANE"]
			samples["RUN"][data["general_informations"][k.split("_")[1]]["RUN"]]["SEQUENCEUR"] = data["general_informations"][k.split("_")[1]]["SEQUENCEUR"]
	nb_sample = len(samples)
	with open(qc, "r") as csv_file:
		qc = {}
		qc["Informations generales"] = {}
		qc["Informations generales"]["Prescription"] = pid_crc
		qc["Informations generales"]["Analyse"] = analyse
		qc["Informations generales"]["Date Analyse"] = date
		qc["Informations generales"]["QC_Version"] = pipe_version
		qc["Informations generales"]["Sample"] = {}
		csv_reader = csv.reader(csv_file, delimiter=',')
		old_row = next(csv_reader)
		qc_order = ["FASTQC STATS", "PICARD COLLECT INSERT SIZE METRICS", "FLAGSTAT STAT", "PICARD MARKEDUPLICATES METRICS", "PICARD Collect Quality Yield Metrics", "PICARD Targeted PCR Metrics"]
		for row in csv_reader:
			for s in samples["ID"]:
				if len(row) > 1 and row[1] == s:
					if not "ID_LABO" in samples["ID"][s]:
						samples["ID"][s]["ID_LABO"] = row[0]
			if len(row) >= 1 and row[0] in qc_order:
				tool = row[0]
				row = next(csv_reader)
				header = row
				qc[tool] = {}
				if tool == "PICARD Targeted PCR Metrics":
					nb_s = 1
				else:
					nb_s = len(samples["ID"])
				for s in range(nb_s):
					row = next(csv_reader)
					stat = row
					lab_id = stat[0]
					id_anon = stat[1]
					qc[tool][samples["ID"][id_anon]['ID_LABO']] = {}
					for h in range(2, len(header)):
						qc[tool][samples["ID"][id_anon]['ID_LABO']][header[h]] = stat[h]
				if tool != "PICARD Targeted PCR Metrics":
					row = next(csv_reader)
			elif len(row) > 0 and row[0] == "PICARD ALIGNMENT METRICS":
				tool = row[0]
				row = next(csv_reader)
				header = row
				qc[tool] = {}
				for k in range(3):
					for s in range(len(samples["ID"])):
						row = next(csv_reader)
						stat = row
						cat = row[2]
						print("CAT:", cat)
						lab_id = stat[0]
						id_anon = stat[1]
						print("id :", samples["ID"][id_anon]['ID_LABO'])
						if not samples["ID"][id_anon]['ID_LABO'] in qc[tool]:
							qc[tool][samples["ID"][id_anon]['ID_LABO']] = {}
						qc[tool][samples["ID"][id_anon]['ID_LABO']][cat] = {}
						for h in range(2, len(header)):
							qc[tool][samples["ID"][id_anon]['ID_LABO']][cat][header[h]] = stat[h]
					row = next(csv_reader)
			elif len(row) > 0 and row[0] in ["GATK DEPTH OF COVERAGE", "GATK CALLABLE LOCI", "PICARD WGS METRICS"]:
				len_sample = 0
				if row[0] == "GATK DEPTH OF COVERAGE" or row[0] == "GATK CALLABLE LOCI":
					flag = False
					for s in samples["ID"]:
						if "WTS" in s:
							flag = True
					if flag == True:
						len_sample = len(samples["ID"]) -1
					if flag == False:
						len_sample = len(samples["ID"])
				if row[0] == "PICARD WGS METRICS":
					flag = False
					for s in samples["ID"]:
						if "WGS-T" in s:
							flag = True
					if flag == True:
						len_sample = 2
					if flag == False:
						len_sample = 1
				if len_sample != 0:
					tool = row[0]
					row = next(csv_reader)
					header = row
					qc[tool] = {}
					for s in range(len_sample):
						row = next(csv_reader)
						stat = row
						lab_id = stat[0]
						id_anon = stat[1]
						qc[tool][samples["ID"][id_anon]['ID_LABO']] = {}
						for h in range(2, len(header)):
							qc[tool][samples["ID"][id_anon]['ID_LABO']][header[h]] = stat[h]
					row = next(csv_reader)
			elif len(row) > 0 and "SNPEFF STAT" in row[0]:
				tool = row[0]
				print(tool)
				header = next(csv_reader)
				qc[tool] = {}
				while "Total/Mean Ratios" not in row:
					row = next(csv_reader)
				stat = row
				print("Stat :", stat[0])
				for s in samples["ID"]:
					qc[tool][samples["ID"][s]["ID_LABO"]] = {}
					for k in range(1, len(header)):
						qc[tool][samples["ID"][s]['ID_LABO']][header[k]] = stat[k]
				row = next(csv_reader)
				ti_tv_header = next(csv_reader)
				ti_tv_data = next(csv_reader)
				for k in range(1, len(ti_tv_header)):
					for s in samples["ID"]:
						if ti_tv_header[k] == s.split("_S")[0]:
							qc[tool][samples["ID"][s]["ID_LABO"]]["Total ratio Ti/Tv"] = ti_tv_data[k]
			elif len(row) > 0 and row[0] == "PICARD VARIANTS METRICS":
				i = []
				for s in samples["ID"]:
					if "WTS" in s:
						i.append(s)
					if "WES-T" in s:
						i.append(s)
					if "WGS-C" in s:
						i.append(s)
				print("i :", i)
				tool = row[0]
				row = next(csv_reader)
				header = row
				row = next(csv_reader)
				stat = row
				qc[tool] = {}
				for s in i:
					qc[tool][samples["ID"][s]["ID_LABO"]] = {}
					for k in range(1, len(header)):
						qc[tool][samples["ID"][s]["ID_LABO"]][header[k]] = stat[k]
				row = next(csv_reader)
				row = next(csv_reader)
				for s in range(len(i)):
					row = next(csv_reader)
					id_anon = row[1]
					hethom = row[2]
					qc[tool][samples["ID"][id_anon]["ID_LABO"]]["Het_hom_ratio"] = hethom
			elif len(row) > 0 and row[0] == "PICARD ILLUMINA BASE CALLING METRICS":
				tool = row[0]
				row = next(csv_reader)
				header = row
				qc["Informations generales"]["PICARD ILLUMINA BASE CALLING METRICS"] = {}
				row = next(csv_reader)
				nb_run = len(samples["RUN"])
				print("NB RUN :", nb_run)
				iteration = 0
				for n in range(nb_run):
					print("L :", row)
					run = row[0]
					nb_lane = len(samples["RUN"][run]["LANE"])
					print("NB_Lane", nb_lane)
					for k in range(nb_lane):
						run = row[0]
						stat = row
						print("S :", stat)
						l = row[1]
						if not run in qc["Informations generales"]["PICARD ILLUMINA BASE CALLING METRICS"]:
							qc["Informations generales"]["PICARD ILLUMINA BASE CALLING METRICS"][run] = {}
						if not l in qc["Informations generales"]["PICARD ILLUMINA BASE CALLING METRICS"][run]:
							qc["Informations generales"]["PICARD ILLUMINA BASE CALLING METRICS"][run][l] = {}
						for k in range(2, len(header)):
							qc["Informations generales"]["PICARD ILLUMINA BASE CALLING METRICS"][run][l][header[k]] = stat[k]
						row = next(csv_reader)
					print("End :", row)
					stat = row
					run = row[1]
					print("Total", run)
					if not "Total" in qc["Informations generales"]["PICARD ILLUMINA BASE CALLING METRICS"][run]:
						qc["Informations generales"]["PICARD ILLUMINA BASE CALLING METRICS"][run]["Total"] = {}
					for k in range(2, len(header)):
						qc["Informations generales"]["PICARD ILLUMINA BASE CALLING METRICS"][run]["Total"][header[k]] = stat[k]
					row = next(csv_reader)
					iteration = iteration + 1
					if iteration != nb_run:
						row = next(csv_reader)
					print("End :", row)
		for s in samples["ID"]:
			print("s :", samples["ID"][s])
			id_labo = samples["ID"][s]["ID_LABO"]
			qc["Informations generales"]["Sample"][id_labo] = {}
			qc["Informations generales"]["Sample"][id_labo]["ID"] = s
			qc["Informations generales"]["Sample"][id_labo]["RUN"] = samples["ID"][s]["RUN"]
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
		json.dump(qc,o, ensure_ascii=False, indent=4, sort_keys = True)
	return output
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Python script which convert an identito csv file into a json')
	parser.add_argument('-q', '--qc_lab', type=str, nargs=1,help='Path of identito csv', required=True)
	parser.add_argument('-p', '--pid_dossier', type=str, nargs=1,help='pid_crc of the prescription', required=True)
	parser.add_argument('-d', '--date_analysis', type=str, nargs=1,help='date of bioinformatic analysis', required=True)
	parser.add_argument('-j', '--pipe_config', type=str, nargs=1,help='path of the pipeline_config', required=True)
	parser.add_argument('-a', '--analysis', type=str, nargs=1,help='analysis id', required=True)
	parser.add_argument('-v', '--version', type=str, nargs=1,help='pipeline version', required=True)
	parser.add_argument('-o', '--output_path', type=str, nargs=1,help='output_file', required=True)
	args = parser.parse_args()
	##
	qc = args.qc_lab[0]
	print("QC :", qc)
	pipe_version = args.version[0]
	pid_crc = args.pid_dossier[0]
	date = args.date_analysis[0]
	json_pipe = args.pipe_config[0]
	print("json :", json_pipe)
	output = args.output_path[0]
	analyse = args.analysis[0]
	##
	qc_json(qc, pid_crc, date, json_pipe, output, analyse, pipe_version)
