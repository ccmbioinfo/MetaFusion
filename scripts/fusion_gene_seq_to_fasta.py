#!/usr/bin/env python
import sys
#sys.path.append("/hpf/largeprojects/ccmbio/jiangyue/DIPG_analysis_by_samples/Scripts/pygeneann/pygenefusionann")
#sys.path.append("/hpf/largeprojects/ccmbio/jiangyue/Genap_ccm/pygenefusionann/")
import pygeneann
import sequtils
import pysam
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('cff_file', action='store', help='CFF file, can be .cff or cff.reann')                                                                                          
parser.add_argument('seq_len', action='store', help='extracted fusion sequence length')

args = parser.parse_args()

seq_len = int(args.seq_len)
for line in open(args.cff_file, "r"):
	fusion = pygeneann.CffFusion(line)
	#annotate fusion id and seq
	print ">"+str(fusion.fusion_id)+"_1"
	print fusion.seq1[-seq_len:]
	print ">"+str(fusion.fusion_id)+"_2"
	print fusion.seq2[:seq_len]

	
