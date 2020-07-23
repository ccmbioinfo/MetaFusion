#!/usr/bin/env python

import sys
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('sum_file', action='store', help='SUMMARY file of fusion reads capture')

args = parser.parse_args()

# SUMMARY GTEX-N7MS-0225-SM-4E3HO Normal DPYD DTD2 NT 0 10 0 20

fusion_stat_dict = {}
read_cnt_thresh = 3
fusion_fraction_thresh = 0.1
for line in open(args.sum_file, "r"):
	tmp = line.split()
	sample = tmp[1]
	disease = tmp[2]
	gene1 = tmp[3]
	gene2 = tmp[4]
	sample_type = tmp[5]
	fusion_cnt = int(tmp[6])
	trans_cnt = int(tmp[7])

	key = (gene1, gene2, disease)

	fusion_stat_dict.setdefault(key, []).append((sample_type, fusion_cnt, trans_cnt))
#print "#FUSION\tTUMOR_CNT\tNORMAL_CNT"
for key in fusion_stat_dict:
	cnts = fusion_stat_dict[key]
	n_tumor = 0
	n_normal = 0
	fusion_fraction = []
	fusion_cnt_list = []
	for sample_type, fusion_cnt, trans_cnt in cnts:
		if fusion_cnt >= read_cnt_thresh:
			if sample_type in ["Tumor", "TP"]:
				n_tumor += 1
			else:
				n_normal += 1
			fusion_fraction.append(float(fusion_cnt)/(fusion_cnt+trans_cnt))
			fusion_cnt_list.append(fusion_cnt)
	if len(fusion_fraction): 
		avg_fusion_fraction = float(sum(fusion_fraction))/len(fusion_fraction)
		max_fusion_fraction = max(fusion_fraction)
	else:
		avg_fusion_fraction = 0
		max_fusion_fraction = 0

	if len(fusion_cnt_list):
		avg_fusion_cnt = float(sum(fusion_cnt_list))/len(fusion_cnt_list)
		max_fusion_cnt = max(fusion_cnt_list)
	else:
		avg_fusion_cnt = 0
		max_fusion_cnt = 0

	print "\t".join(map(str, [" ".join(list(key)), n_tumor, n_normal, avg_fusion_cnt, max_fusion_cnt, avg_fusion_fraction, max_fusion_fraction, len(cnts)]))

	
