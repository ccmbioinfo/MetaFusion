#!/usr/bin/env python
import sys
import pygeneann
import sequtils
import pysam
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('cff_file', action='store', help='CFF file, can be .cff or cff.reann')
parser.add_argument('ensbed', action='store', help='Ensemble gene file')
parser.add_argument('ref_fa', action='store', help='Reference genome file')

args = parser.parse_args()

gene_ann = pygeneann.GeneAnnotation(args.ensbed)
	
ref = pysam.FastaFile(args.ref_fa)
seq_dict ={}
for line in open(args.cff_file, "r"):
	fusion = pygeneann.CffFusion(line)

	gene1 = fusion.reann_gene1
	gene2 = fusion.reann_gene2
	fusion_id = fusion.fusion_id

	trans_seqs1 = pygeneann.build_transcript_sequences(gene_ann, fusion.reann_gene1, ref)	
	trans_seqs2 = pygeneann.build_transcript_sequences(gene_ann, fusion.reann_gene2, ref)	

	if not trans_seqs1 or not trans_seqs2:
		continue

	for trans_id in trans_seqs1:
		trans_seq = trans_seqs1[trans_id]
		print ">Gene1_" + gene1 + "_" + gene2 + "_" + trans_id + "_" + fusion_id
		print trans_seq
	for trans_id in trans_seqs2:
		trans_seq = trans_seqs2[trans_id]
		print ">Gene2_" + gene1 + "_" +  gene2 + "_" + trans_id + "_" + fusion_id
		print trans_seq
#fasta_f.close()
