#!/usr/bin/env python

import sys
#sys.path.append("/hpf/largeprojects/ccmbio/jiangyue/DIPG_analysis_by_samples/Scripts/pygeneann/pygenefusionann")
import pysam
import sequtils
import os
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('bam_file', action='store', help='BAM file from which to extract fusion-supporting pairs')

args = parser.parse_args()


bam_file = sys.argv[1]

bam = pysam.AlignmentFile(args.bam_file, "rb")
ref = bam.references[0]
tmp = ref.split("_")
bp = int(tmp[-1]) - 1

def get_softclip_len(read):
	sc_len_head = 0
	sc_len_tail = 0
	t = read.cigartuples[0]
	if t[0] == 4: # softclipping
		sc_len_head += int(t[1])
	t = read.cigartuples[-1]
	if t[0] == 4: # softclipping
		sc_len_tail += int(t[1])
	return sc_len_head, sc_len_tail, sc_len_head+sc_len_tail
filename = os.path.basename(args.bam_file)
dirname = os.path.dirname(args.bam_file)
tmp = filename.split(".")
sample = tmp[0]
outfa = open(args.bam_file + ".fa", "w")
for refname in bam.references:
	captured_reads_dict = {}
	tmp = refname.split("_")
	gene1 = tmp[0]
	gene2 = tmp[1]
	outseqs = []
	for read in bam.fetch(refname):
		if read.is_unmapped or read.is_secondary or read.is_supplementary:
			continue
		if read.reference_start	< bp -10 and read.reference_end >= bp and get_softclip_len(read)[2] <= read.query_length * 0.1: # potential fusion/transcript support reads, span >=10 bp junction sequences, and <= 10% clipping
			read_id = "1" if read.is_read1 else "2"
			print >> outfa, ">" + read.query_name + "_" + read_id
			print >> outfa, read.query_sequence
			#captured_reads_dict.setdefault(read.query_name, "")
		

outfa.close()
