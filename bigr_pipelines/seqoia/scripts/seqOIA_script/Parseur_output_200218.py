#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Mario NEOU
# 2020/02/18
#python3


# /home/mneou/Documents/Projets/Scripts/Arriba2FusionInspector/MAP567_WTS_fusions.tsv


import re
import sys



dict_file=sys.argv[1]
print(dict_file)

output = sys.argv[2]
print(output)
#dict_file="/home/mneou/Documents/Projets/Scripts/Arriba2FusionInspector/MAP567_WTS_fusions.tsv"

data = []
with open(dict_file, 'r') as f:
	f_list = f.readlines()
	for lines in f_list:
		lines = lines.split("\t")
		partenaires1 = lines[0]
		partenaires2 = lines[1]
		print("partenaires :", partenaires1, partenaires2)
		list_partenaire1 = partenaires1.split(",")
		list_partenaire2 = partenaires2.split(",")
		print("list_partenaire2 :", list_partenaire2)
		print("list_partenaire1 :", list_partenaire1)
		for i in range(len(list_partenaire2)):
			list_partenaire2[i] = re.sub('\([0-9]+\)','',list_partenaire2[i])
		for i in range(len(list_partenaire1)):
			list_partenaire1[i] = re.sub('\([0-9]+\)','',list_partenaire1[i])
		print("list_partenaire1_modif :", list_partenaire1)
		print("list_partenaire2_modif :", list_partenaire2)
		for p1 in list_partenaire1:
			for p2 in list_partenaire2:
				data.append([p1,p2])
				print("FINAL :", [partenaires1, partenaires2])

#dict_file.split(".tsv")
#output = dict_file.split(".tsv")[0] +"_FusionInspector.tsv"

f = open(output, "w")
for u in data :
    #f.write( {}\t{} )format(u[0],u[1])
	f.write( '--'.join(u) )
	f.write( '\n' )
f.close()


            
        
		#if lines[0] == "@SQ":
			#chromosome = lines[1].split("SN:")[1]
			#end = lines[2].split("LN:")[1]
			#data.append([chromosome, end])
                
            
"""
def dict_data_extraction(dict_file):
    dict_file="/home/mneou/Documents/Projets/Scripts/Arriba2FusionInspector/MAP567_WTS_fusions.tsv"
	data = []
	with open(dict_file, 'r') as f:
		f_list = f.readlines()
		for lines in f_list:
			lines = lines.split("\t")
			if lines[0] == "@SQ":
				chromosome = lines[1].split("SN:")[1]
				end = lines[2].split("LN:")[1]
				data.append([chromosome, end])
	return data
def bed_generation(data,output):
	with open(output, "w") as f:
		for k in data:
			line =[str(k[0]), "0", str(k[-1]) + "\n"]
			line = "\t".join(line)
			f.write(line)
	print(output)
	return output

if __name__ == '__main__':
	dict_file = "/home/alegendre/Script/GRCh38.92.dict"
	output = "/home/alegendre/Script/GRCH38.92.bed"
	data = dict_data_extraction(dict_file)
	output = bed_generation(data,output)

"""
