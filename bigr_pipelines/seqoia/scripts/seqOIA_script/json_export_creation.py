# -*- coding: utf-8 -*-
# Adrien Legendre 
# 2019/08/26
from __future__ import unicode_literals
import json
import os
import argparse
import sys
import glob

# this script will add the path of the vcf, cluster_config and pipeline_config into the json
# Purity, tmb, msi, infered ploidy extraction
def genome_infos(facet, msi_f):
	with open(facet, "r") as f:
		lines = f.readline()
		lines = next(f)
		purity = lines.split("\t")[15]
		ploidy = lines.split("\t")[16]
		ploidy = round(float(ploidy),1)
	with open(msi_f, "r") as f:
		lines = f.readline()
		lines = next(f)
		msi = lines.split("\t")[-1].split("\n")[0]
	return purity, ploidy, msi

def tmb_infos(tmb_f):
	with open(tmb_f, "r") as f:
		for lines in f:
			if "TMB=" in lines:
				tmb = lines.split("TMB= ")[-1].split("\n")[0]
	return tmb

#extraction of pipeline version form pipeline_config:
def pipeline_version_extraction(config):
	config_file = glob.glob(config + "/*")
	if len(config_file) != 2:
		print("Configuration files (cluster and pipeline) are not present in irods")
		return sys.exit(1)
	for k in config_file:
		if "cluster" in k:
			cluster_config = k
		if "pipeline" in k :
			pipeline_config = k
	with open(pipeline_config, "r", encoding='utf-8') as f:
		data = json.load(f)
		id_anon_extended = data["general_informations"]["SAMPLES"]
		pipeline_version = data["general_informations"]["ID_PIPELINE"]
		version = data["general_informations"]["VCF_VERSION"]
		compte_rendu = data["Compte_rendu"]
		if data["general_informations"]["PIPELINE_CONFIG"] == "PIPELINE_CONFIG_VARIABLE" or data["general_informations"]["CLUSTER_CONFIG"] == "CLUSTER_CONFIG_VARIABLE":
			return "error"
		cluster_config_path = os.path.abspath(cluster_config)
		pipeline_config_path = os.path.abspath(pipeline_config)
		general_info = data["general_informations"]
	return pipeline_version, cluster_config_path, pipeline_config_path, id_anon_extended, compte_rendu, version, general_info

#lecture of the json
def json_creation(json_file,info_template, info_variable,id_anon_extended, output, compte_rendu, version, general_info, bam,bai):
	with open(json_file, "r", encoding='utf-8') as f:
		data = json.load(f)
	if len(info_template) == len(info_variable):
		for x,y in zip(info_template, info_variable):
			data[x] = y
		data["run_id"] = {}
		data["bam"] = {}
		for k in id_anon_extended:
			if "WGS-C" in k:
				data["WGS-C_extended"] = k.split("_")[2]
				data["run_id"]["WGS-C"] = general_info["WGS-C"]["RUN"]
			if "WES-T" in k:
				data["WES-T_extended"] = k.split("_")[2]
				data["run_id"]["WES-T"] = general_info["WES-T"]["RUN"]
			if "WTS" in k:
				data["WTS_extended"] = k.split("_")[2]
				data["run_id"]["WTS"] = general_info["WTS"]["RUN"]
			if "WGS-T" in general_info:
				data["run_id"]["WGS-T"] = general_info["WGS-T"]["RUN"]
			data["bam"][k.split("_S")[0]] = {}
			for chrom in ["1","2","3","4","5","6","7","8","9","10","11","12","13", "14","15","16","17","18","19","20","21","22","X","Y","MT"]:
				data["bam"][k.split("_S")[0]][chrom] = {}
				for b in bam:
					if "WTS" in k:
						if "{ids}_chr_{chrom}_wts_markdup.bam".format(ids=k, chrom=chrom) in b:
							data["bam"][k.split("_S")[0]][chrom]["bam"] = "fs://" + b
					else:
						if "{ids}_chr_{chrom}_markdup.bam".format(ids=k, chrom=chrom) in b:
							data["bam"][k.split("_S")[0]][chrom]["bam"] = "fs://" + b
				for i in bai:
					if "WTS" in k:
						if "{ids}_chr_{chrom}_wts_markdup.bam.bai".format(ids=k, chrom=chrom) in i:
							data["bam"][k.split("_S")[0]][chrom]["bai"] = "fs://" + i
					else:
						if "{ids}_chr_{chrom}_markdup.bam.bai".format(ids=k, chrom=chrom) in i:
							data["bam"][k.split("_S")[0]][chrom]["bai"] = "fs://" + i
		data["TAR"] = info_variable[5]
		data["cnv_fusion"] = info_variable[5]
		data["CNV"] = info_variable[6]
		data["vcf_version"] = version
		data["Compte_rendu"] = compte_rendu
		with open(output, "w", encoding = 'utf-8') as f:
			json.dump(data,f, ensure_ascii=False)
		return output
	else:
		return "error in list, some informations are missing (vcf, pipeline configuration, cluster configuration or pipeline version)"


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Python script which extract all informations about tools versions and options, pipeline version and vcf path in order to create a json')
	parser.add_argument('-v', '--vcf', type=str, nargs=1,help='Path of vcf', required=True)
	parser.add_argument('-j', '--json', type=str, nargs=1,help='Path of json file for the vcf for the proband', required=True)
	parser.add_argument('-o', '--output',type=str, nargs=1,help='Path of the output file', required=True)
	parser.add_argument('-q', '--qc',type=str, nargs=1,help='Path of the qc file', required=True)
	parser.add_argument('-c', '--config',type=str, nargs=1,help='Path of the config folder in irods', required=True)
	parser.add_argument('-cnv', '--cnv_path',type=str, nargs=1,help='Path of the config cnv folder in irods', required=False)
	parser.add_argument('-tar', '--tar_path',type=str, nargs=1,help='Path of the config cnv folder in irods', required=False)
	parser.add_argument('-tmb', '--tmb_file',type=str, nargs=1,help='Path of tmb file', required=False)
	parser.add_argument('-msi', '--msi_file',type=str, nargs=1,help='Path of msi file', required=False)
	parser.add_argument('-f', '--facet_bed',type=str, nargs=1,help='Path of facet bed', required=False)
	parser.add_argument('-bam', '--bam_path',type=str, nargs='+',help='Path of the bam files', required=True)
	parser.add_argument('-bai', '--bai_path',type=str, nargs='+',help='Path of the bai files', required=True)
#########################
	args = parser.parse_args()
	vcf = args.vcf[0]
	config = args.config[0]
	output = args.output[0]
	json_file = args.json[0]
	qc = args.qc[0]
	bam = args.bam_path
	bai = args.bai_path
	if args.tmb_file:
		tmb_f = args.tmb_file[0]
	else:
		tmb = "Not Present"
	if args.msi_file:
		msi_f = args.msi_file[0]
	else:
		msi = "Not Present"
	if args.facet_bed:
		facet = args.facet_bed[0]
	else:
		purity = "No Present"
		ploidy = "No Present"
	if args.tar_path:
		tar = args.tar_path[0]
	else:
		tar = "Not present"
	if args.cnv_path:
		cnv = args.cnv_path[0]
	else:
		cnv = "Not Present"
	########################
	if pipeline_version_extraction(config) == "error":
		print("PIPELINE_CONFIG and CLUSTER_CONFIG_VARIABLE are still default variable, please manually change them")
		sys.exit(1)
	else:
		pipe_version = pipeline_version_extraction(config)[0]
		cluster_config_path = pipeline_version_extraction(config)[1]
		pipeline_config_path = pipeline_version_extraction(config)[2]
		id_anon_extended = pipeline_version_extraction(config)[3]
		compte_rendu = pipeline_version_extraction(config)[4]
		version = pipeline_version_extraction(config)[5]
		if args.facet_bed:
			purity  = genome_infos(facet, msi_f)[0]
			ploidy = genome_infos(facet, msi_f)[1]
		if args.msi_file:
			msi = genome_infos(facet, msi_f)[2]
		if args.tmb_file:
			tmb = tmb_infos(tmb_f)
		general_info = pipeline_version_extraction(config)[6]
		#building of the python dictionnary in order to introduce it to json
		info_template = ['vcf', 'pipeline configuration', 'cluster configuration', 'pipeline version', 'qc', 'tar_output', 'cnv', 'purity', 'ploidy', 'tmb', 'msi']
		info_variable = [vcf, pipeline_config_path, cluster_config_path, pipe_version, qc, tar, cnv, purity, str(ploidy), tmb, msi]
		#Write the dictionnary into the json
		json_output = json_creation(json_file,info_template, info_variable,id_anon_extended, output, compte_rendu, version, general_info,bam,bai)
		print(json_output)
