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


bam = pysam.AlignmentFile(args.bam_file, "rb")
ref = bam.references[0]
tmp = ref.split("_")
#filename = (args.bam_file.split("/")[-1])
filename = os.path.basename(args.bam_file)
dirname = os.path.dirname(args.bam_file)
tmp = filename.split(".")
sample = tmp[0]
outfa = open(args.bam_file + ".fa", "w")

pair_read_dict = {}
for read in bam.fetch(until_eof=True):
	if read.is_unmapped or read.is_secondary or read.is_supplementary or read.mate_is_unmapped:
		continue
	tmp1 = read.reference_name.split("_")
	tmp2 = read.next_reference_name.split("_")

	gene_id = tmp1[0]
	mate_gene_id = tmp2[0]
	fusion_id = tmp1[-1]
	mate_fusion_id = tmp2[-1]
	
	nm = read.get_tag("NM")
	if gene_id != mate_gene_id and fusion_id == mate_fusion_id and nm <= 5:
		pair_read_dict.setdefault(read.query_name, []).append(read)	
pairedreads = pysam.AlignmentFile(args.bam_file + ".filtered.sam", "w", template=bam)
for rname in pair_read_dict:
	reads = pair_read_dict[rname]
	if len(reads) == 2:
		pairedreads.write(reads[0])
		pairedreads.write(reads[1])
		
		read = reads[0]
		read_id = "1" if read.is_read1 else "2"
		print >> outfa, ">" + read.query_name + "_" + read_id
		print >> outfa, read.query_sequence

		read = reads[1]
		read_id = "1" if read.is_read1 else "2"
		print >> outfa, ">" + read.query_name + "_" + read_id
		print >> outfa, read.query_sequence
pairedreads.close()	
outfa.close()
