#!/usr/bin/env python2
# -*- coding: utf-8 -*-
# Adrien Legendre 
# 2020/09/22
#python2
# Modificator of CNV result from WisecondorX and Facet and fusion

import argparse
import glob
import json
import sys
import os

def json_creator(mercury, census, curie):
	json_list = {}
	with open(mercury, "r") as f:
		data = f.readlines()
		for lines in data:
			lines = lines.split("\r")[0].split("\n")
			if lines[0] != '':
				json_list[lines[0]] = ["Mercury","1","NA","0","NA","0","NA"]
	with open(census, "r") as f:
		data = f.readlines()
		for lines in data:
			lines = lines.split(",")
			if lines[0] in json_list:
				json_list[lines[0]][4] = "Census"
				json_list[lines[0]][5] = "1"
				json_list[lines[0]][6] = lines[1] + "|" + lines[2]
			else:
				json_list[lines[0]] = ["NA", "0", "NA", "0", "Census","1", lines[1] + "|" + lines[2]]
	with open(curie, "r") as f:
		data = f.readlines()
		for lines in data:
			lines = lines.split("\n")[0]
			if lines in json_list:
				json_list[lines][2] = "Curie"
				json_list[lines][3] = "1"
			else:
				json_list[lines] = ["NA","0", "Curie","1", "NA","0", "NA"]
	return json_list

def CNV_modificator(facet_genome, facet_chr, wisecondor, json_list):
	facet = [facet_genome] + facet_chr
	for k in facet:
		data_facet = []
		print("FICHIER INPUT", k)
		with open(k, "r") as f:
			for lines in f:
				l = lines.split("\t")
				if l[0] == "AnnotSV ID":
					output_l = l[:21] + ["Mercury_list", "%Mercury_list", "Curie_list", "%Curie_list","Gene_Census_list", "%Gene_Census_list", "%All_list", "Tier|Role"] + l[21:49] + l[50:60] + l[61:]
					output_l[4] = "SV " + output_l[4]
					output_l[5] = output_l[5] + " (Mb)"
				else:
					l[5] = str(float(l[5])/1000000)
					if l[19] == "split":
						if str(l[20]) in json_list:
							all_list = ["0"]
							if "1" in json_list[str(l[20])]:
								all_list = ["1"]
							output_l = l[:21] + json_list[str(l[20])][:-1] + all_list + [json_list[str(l[20])][-1]] + l[21:49] + l[50:60] + l[61:]
						else:
							output_l = l[:21] + ["NA","0", "NA","0", "NA","0", "0", "NA"] + l[21:49] + l[50:60] + l[61:]
					elif l[19] == "full":
						genes = l[20].split("/")
						genes = list(set(genes))
						len_g = float(len(genes))
						all_list = 0
						tr = []
						p_mercury = 0
						p_curie = 0
						p_census = 0
						for g in genes:
							if g in json_list:
								all_list = all_list + 1
								if "Mercury" in json_list[g]:
									p_mercury = p_mercury + 1
								if "Curie" in json_list[g]:
									p_curie = p_curie + 1
								if "Census" in json_list[g]:
									p_census = p_census + 1
									tr.append(json_list[g][-1])
								else:
									tr.append(json_list[g][-1])
							if g not in json_list:
								tr.append("NA")
						output_l = l[:21] + [str(p_mercury),str((float(p_mercury)/len_g)*100),str(p_curie),str((float(p_curie)/len_g)*100),str(p_census),str((float(p_census)/len_g)*100), str((float(all_list)/len_g)*100), "/".join(tr)] + l[21:49] + l[50:60] + l[61:]
					else:
						output_l = l[:21] + ["NA","0","NA","0","NA", "0", "0","NA"] + l[21:49] + l[50:60] + l[61:]
				data_facet.append(output_l)
		output_facet = k.replace("annotated.tsv", "annotated.final.tsv")
		with open(output_facet, "w+") as o:
			print("OUTPUT :", output_facet)
			for k in data_facet:
				k = k[:4] + k[5:]
				o.write("\t".join(k))
	with open(wisecondor, 'r') as f:
		data_wise = []
		print("FICHIER INPUT", wisecondor)
		for lines in f:
			l = lines.split("\t")
			if l[0] == "AnnotSV ID":
				output_l = l[:5]+ ["Ratio", "Zscore", "Type"] + l[8:10] + ["Mercury_list", "%Mercury_list", "Curie_list", "%Curie_list","Gene_Census_list", "%Gene_Census_list", "%All_list", "Tier|Role"] + l[10:38] + l[39:49] + l[50:]
			else:
				if l[8] == "split":
					if str(l[9]) in json_list:
						all_list = ["0"]
						if "1" in json_list[str(l[9])]:
							all_list = ["1"]
						output_l = l[:10] + json_list[str(l[9])][:-1] + all_list + [json_list[str(l[9])][-1]] + l[10:38] + l[39:49] + l[50:]
					else:
						output_l = l[:10] + ["NA", "0", "NA", "0", "NA", "0", "0", "NA"] + l[10:38] + l[39:49] + l[50:]
				elif l[8] == "full":
					genes = l[9].split("/")
					genes = list(set(genes))
					len_g = float(len(genes))
					all_list = 0
					tr = []
					p_mercury = 0
					p_curie = 0
					p_census = 0
					for g in genes:
						if g in json_list:
							all_list = all_list + 1
							if "Mercury" in json_list[g]:
								p_mercury = p_mercury + 1
							if "Curie" in json_list[g]:
								p_curie = p_curie + 1
							if "Census" in json_list[g]:
								p_census = p_census + 1
								tr.append(json_list[g][-1])
							else:
								tr.append(json_list[g][-1])
						if g not in json_list:
							tr.append("NA")
					output_l = l[:10] + [str(p_mercury),str((float(p_mercury)/len_g)*100),str(p_curie),str((float(p_curie)/len_g)*100),str(p_census),str((float(p_census)/len_g)*100), str((float(all_list)/len_g)*100), "/".join(tr)] + l[10:38] + l[39:49] + l[50:]
				else:
					output_l = l[:10] + ["NA","0","NA","0","NA", "0", "0","NA"] + l[10:38] + l[39:49] + l[50:]
				if int(l[3]) >= int(l[2]):
					output_l[4] = str((float(l[3]) - float(l[2]))/1000000) + " Mb"
				else:
					output_l[4] = str((float(l[2]) - float(l[3]))/1000000) + " Mb"
			data_wise.append(output_l)
		output_wisecondor = wisecondor.replace("annotated.tsv", "annotated.final.tsv")
		with open(output_wisecondor, "w+") as o:
			print("OUTPUT :", output_wisecondor)
			for k in data_wise:
				o.write("\t".join(k))
	return wisecondor

def json_fusion(arriba, star, catcher):
	data = {}
	with open(arriba, "r") as f:
		for lines in f:
			lines = lines.split("\t")
			ids = [str(lines[0]),  str(lines[1])]
			if "#gene1" in lines[0]:
				continue
			else:
				for k in  range(len(ids)):
					if  ","  in  ids[k]:
						ids[k]  =  ids[k].split(",")
						for  i  in  range(len(ids[k])):
							ids[k][i]  =  ids[k][i].split("(")[0]
					else:
						if  "("  in  ids[k]:
							ids[k]  =  ids[k].split("(")[0]
				if  type(ids[0])  ==  str  and  type(ids[1])  ==  str:
					if ids[0] + "--" + ids[1] not in data:
						data[ids[0] + "--" + ids[1]] = [[str(int(lines[11]) + int(lines[12]))],[lines[13]], "NA", "NA", "NA", "NA"]
					else:
						data[ids[0] + "--" + ids[1]][1].append(lines[13])
						data[ids[0] + "--" + ids[1]][0].append(str(int(lines[11]) + int(lines[12])))
				elif  type(ids[0])  ==  str  and  type(ids[1])  ==  list:
					for  k  in  ids[1]:
						if ids[0] + "--" + k not in data:
							data[ids[0] + "--" + k] = [[str(int(lines[11]) + int(lines[12]))], [lines[13]], "NA", "NA", "NA", "NA"]
						else:
							data[ids[0] + "--" + k][1].append(lines[13])
							data[ids[0] + "--" + k][0].append(str(int(lines[11]) + int(lines[12])))
				elif  type(ids[0])  ==  list  and  type(ids[1])  ==  str:
					for  k  in  ids[0]:
						if k + "--" + ids[1] not in data:
							data[k + "--" + ids[1]] = [[str(int(lines[11]) + int(lines[12]))], [lines[13]], "NA", "NA", "NA", "NA"]
						else:
							data[k + "--" + ids[1]][1].append(lines[13])
							data[k + "--" + ids[1]][0].append(str(int(lines[11]) + int(lines[12])))
				else:
					for  k  in  ids[0]:
						for  x  in  ids[1]:
							if k + "--" + x  not in data:
								data[k + "--" + x] = [[str(int(lines[11]) + int(lines[12]))],[lines[13]], "NA", "NA", "NA", "NA"]
							else:
								data[k + "--" + x][1].append(lines[13])
								data[k + "--" + x][0].append(str(int(lines[11]) + int(lines[12])))
	with open(star, "r") as f:
		for lines in f:
			lines = lines.split("\t")
			if "#FusionName" in lines[0]:
				continue
			else:
				if lines[0] not in data:
					data[lines[0]] = ["NA","NA",[lines[1]], [lines[2]], "NA", "NA"]
				else:
					if data[lines[0]][2] == "NA":
						data[lines[0]][2] =[lines[1]]
					else:
						data[lines[0]][2].append(lines[1])
					if data[lines[0]][3] == "NA":
						data[lines[0]][3] = [lines[2]]
					else:
						data[lines[0]][3].append(lines[2])
	with open(catcher, "r") as f:
		for lines in f:
			lines = lines.split("\t")
			ids = lines[0] + "--" + lines[1]
			if "Gene_1_symbol" in lines[0]:
				continue
			else:
				if ids not in data:
					data[ids] = ["NA","NA","NA", "NA", [lines[4]], [lines[5]]]
				else:
					if data[ids][4] == "NA":
						data[ids][4] = [lines[4]]
					else:
						data[ids][4].append(lines[4])
					if data[ids][5] == "NA":
						data[ids][5] = [lines[5]]
					else:
						data[ids][5].append(lines[5])
	return data

def Fusion_modificator(fusion, json_fusion, json_CNV):
	data = []
	with open(fusion, 'r') as f:
		print(fusion)
		for lines in f:
			lines = lines.split("\t")
			if "#FusionName" in lines[0]:
				output_l = lines[0:5] + ["Arriba", "Star-fusion", "Fusion-catcher", "Fusion_tools_count","Mercury_Left_Gene", "Curie_Left_Gene", "Census_Left_Gene", "GeneCensus_Tier_Role_Left_Gene", "Mercury_Right_Gene", "Curie_Right_Gene", "Census_Right_Gene", "GeneCensus_Tier_Role_Right_Gene", "%Fusion_Mercury_list", "%Fusion_Curie_list", "%Fusion_Census_list", "%Fusion_All_list", "Census_Tier_Role_fusion", "Arriba_Junction_Reads",  "Arriba_Spanning_Read",  "Star_fusion_Junction_Reads",  "Star_Fusion_Spanning_Read",  "Fusion_Catcher_Spanning_pairs",  "Fusion_Catcher_Spanning_unique_reads"]  +  lines[8:]
			else:
				right = lines[0].split("--")[1]
				left = lines[0].split("--")[0]
				p_mercury = 0
				p_curie = 0
				p_census = 0
				p_all = 0
				tr = ["NA", "NA"]
				if right in json_CNV:
					if "Mercury" in json_CNV[right]:
						p_mercury = p_mercury + 0.5
					if "Curie" in json_CNV[right]:
						p_curie = p_curie + 0.5
					if "Census" in json_CNV[right]:
						p_census = p_census + 0.5
				if left in json_CNV:
					if "Mercury" in json_CNV[left]:
						p_mercury = p_mercury + 0.5
					if "Curie" in json_CNV[left]:
						p_curie = p_curie + 0.5
					if "Census" in json_CNV[left]:
						p_census = p_census + 0.5
				if right not in json_CNV and left in json_CNV:
					p_all = 0.5
					tr[0] = json_CNV[left][-1]
					infos_cnv = [json_CNV[left][0],json_CNV[left][2],json_CNV[left][4],json_CNV[left][-1], "NA", "NA", "NA", "NA"]
				elif right in json_CNV and left not in json_CNV:
					p_all = 0.5
					tr[1] = json_CNV[right][-1]
					infos_cnv = ["NA", "NA", "NA", "NA", json_CNV[right][0],json_CNV[right][2],json_CNV[right][4], json_CNV[right][-1]]
				elif right in json_CNV and left in json_CNV:
					p_all = 1
					tr = [json_CNV[left][-1], json_CNV[right][-1]]
					infos_cnv = [json_CNV[left][0],json_CNV[left][2],json_CNV[left][4],json_CNV[left][-1], json_CNV[right][0],json_CNV[right][2],json_CNV[right][4],json_CNV[right][-1]]
				else:
					infos_cnv = ["NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA"]
				new_infos = [str(p_mercury), str(p_curie), str(p_census), str(p_all),"/".join(tr)]
				if lines[0] in json_fusion:
					infos_fusion = json_fusion[lines[0]]
					for k in range(len(infos_fusion)):
						if type(infos_fusion[k]) == list and len(infos_fusion[k]) == 1:
							infos_fusion[k] = str(infos_fusion[k][0])
						elif type(infos_fusion[k]) == list and len(infos_fusion[k]) > 1:
							infos_fusion[k] = "/".join(infos_fusion[k])
						else:
							continue
					tools = ["Absent", "Absent", "Absent", 0]
					if infos_fusion[0] != "NA"  and infos_fusion[1] != "NA":
						tools[0] = "Present"
						tools[-1] = tools[-1] + 1
					if infos_fusion[2] != "NA"  and infos_fusion[3] != "NA":
						tools[1] = "Present"
						tools[-1] = tools[-1] + 1
					if infos_fusion[4] != "NA"  and infos_fusion[5] != "NA":
						tools[2] = "Present"
						tools[-1] = tools[-1] + 1
					tools[-1] = str(tools[-1])
				else:
					tools = ["Absent", "Absent", "Absent", "In silico reciprocal"]
					infos_fusion = ["NA", "NA", "NA", "NA", "NA", "NA"]
				output_l = lines[0:5] + tools + infos_cnv + new_infos + infos_fusion + lines[8:]
			data.append(output_l)
	output_fusion = fusion.replace("tsv.coding_effect", "coding_effect.final.tsv")
	with open(output_fusion, "w+") as o:
		for k in range(len(data)):
			for x in data[k]:
				if type(x) == int:
					print("ERROR:",x)
			o.write("\t".join(data[k]))
	return "Done"

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='VDM transcript analyse')
	parser.add_argument('-m', '--mercury_list', type=str, nargs=1,action='append', default = ["/scratch2/tmp/shared_files_tmp/CNV/gene_list/cancer_gene_mercury.list.txt"], help='Path of output mercury gene list', required=False)
	parser.add_argument('-c', '--curie_list', type=str, nargs=1,action='append', default = ["/scratch2/tmp/shared_files_tmp/CNV/gene_list/05-06-2018_WES_CancerGenesList_Large.txt"], help='Path of output curie gene list', required=False)
	parser.add_argument('-g', '--gene_census_list', type=str, nargs=1,action='append', default = ["/scratch2/tmp/shared_files_tmp/CNV/gene_list/cancer_gene_census_list_annotated.txt"], help='Path of output census gene list', required=False)
	parser.add_argument('-fc', '--facet_chr_file', type=str, nargs="+",action='append', help='Path of files from AnnotSV facet chromosome', required=False)
	parser.add_argument('-fg', '--facet_genome_file', type=str, nargs=1,action='append', help='Path of file from AnnotSV facet genome', required=False)
	parser.add_argument('-w', '--wisecondor_file', type=str, nargs=1,action='append', help='Path of file from AnnotSV Wisecondor genome', required=False)
	parser.add_argument('-fi', '--fusion_inspector_file', type=str, nargs=1,action='append', help='Path of file from Fusion Inspector', required=False)
	parser.add_argument('-a', '--arriba_file', type=str, nargs=1,action='append', help='Path of file from Arriba', required=False)
	parser.add_argument('-f', '--fusion_catcher_file', type=str, nargs=1,action='append', help='Path of file from Fusion Catcher', required=False)
	parser.add_argument('-s', '--star_file', type=str, nargs=1,action='append', help='Path of file from Star', required=False)
	###########
	args = parser.parse_args()
	mercury = args.mercury_list[0]
	census = args.gene_census_list[0]
	curie = args.curie_list[0]
	if args.facet_genome_file:
		facet_genome = args.facet_genome_file[0][0]
		facet_chr = args.facet_chr_file[0]
		wisecondor = args.wisecondor_file[0][0]
	if args.fusion_inspector_file:
		fusion = args.fusion_inspector_file[0]
		arriba = args.arriba_file[0]
		catcher = args.fusion_catcher_file[0]
		star = args.star_file[0]
	#####
	if args.fusion_inspector_file:
		print("Arriba:", arriba)
		print("Star: ", star)
		print("Catcher :", catcher)
		json_fusion = json_fusion(arriba[0], star[0], catcher[0])
		json_CNV = json_creator(mercury, census, curie)
		Fusion_modificator(fusion[0], json_fusion, json_CNV)
	if args.facet_genome_file:
		print("CNV")
		print("FACET GENOME:", facet_genome)
		print("FACET CHR:", facet_chr)
		json_CNV = json_creator(mercury, census, curie)
		print("JSON FINISH")
		CNV_modificator(facet_genome, facet_chr, wisecondor, json_CNV)
		print("CNV MODIF")
