#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Adrien Legendre
# SeqOIA
# 21/06/19
# Allow to sum up the identito-vigilance file for SeqOIA runs

import glob
import csv
import os
import argparse
import datetime

parser = argparse.ArgumentParser(description='Create the identito-vigilance resume for a SeqOIA run with all the samples')
parser.add_argument('-i', '--sample_identito_csv', type=str, nargs='*', action='append', help='csv file from identito-vigilance.py for all samples', required=True)
parser.add_argument('-o', '--output', type=str, nargs=1,action='store', help='csv output file', required=True)
args = parser.parse_args()

identito = args.sample_identito_csv[0]
output = args.output[0]


pwd = os.getcwd()

data =[]
for k in identito:
	with open(k, "r") as csv_file:
		list_file = []
		csv_file = csv.reader(csv_file)
		for row in csv_file:
			list_file.append(row)
		data.append(list_file[-1])

with open(output, "w+") as csv_file:
	writer = csv.writer(csv_file, quoting=csv.QUOTE_ALL)
	writer.writerow(["Run","Date"])
	writer.writerow([output.split('/')[-2], datetime.datetime.now()])
	writer.writerow("")
	writer.writerow(["Samples", "Number of Concordant", "Number of Inconcordant", "Number of no call", "Homologie"])
	for k in data:
		writer.writerow(k)
