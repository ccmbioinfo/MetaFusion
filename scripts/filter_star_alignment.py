#!/usr/bin/env python

import sys
import os
import commands
import argparse
import pysam

#sys.path.append("/hpf/largeprojects/ccmbio/jiangyue/DIPG_analysis_by_samples/Scripts/pygeneann/pygenefusionann")
import pygeneann
import sequtils

class Cluster():
	def __init__(self, chr1, pos1, chr2, pos2, n):
		self.chr1 = chr1
		self.pos1 = pos1
		self.chr2 = chr2
		self.pos2 = pos2
		self.number = n
	def tostring(self):
		return "\t".join(map(lambda x:str(x), [self.chr1, self.pos1, self.chr2, self.pos2, self.number]))

def load_bam_dict(bam_list):
	bam_dict = {}
	for line in open(bam_list, "r"):
		tmp = line.split()
		sample = tmp[0]
		bamfile = tmp[1]
		bam_dict.setdefault(sample, bamfile)
	return bam_dict

def format_chr(ref):
	if not ref.startswith("chr"):
		return "chr" + ref
	else:
		return ref
	

def output_paired_fastq(read1, read2, star_bam, fastq1, fastq2):
	sam_line = read1.tostring(star_bam).strip()
	base_qualities = sam_line.split()[10]
	print >> fastq1, "@"+read1.query_name+"/1"
	if read1.is_reverse:
		print >> fastq1, sequtils.rc_seq(read1.query_sequence, "rc")
		print >> fastq1, "+"
		print >> fastq1, sequtils.rc_seq(base_qualities, "r")		
	else:
		print >> fastq1, read1.query_sequence
		print >> fastq1, "+"
		print >> fastq1, base_qualities

	sam_line = read2.tostring(star_bam).strip()
	base_qualities = sam_line.split()[10]
	print >> fastq2, "@"+read1.query_name+"/2"
	
	if read2.is_reverse:
		print >> fastq2, sequtils.rc_seq(read2.query_sequence, "rc")
		print >> fastq2, "+"
		print >> fastq2, sequtils.rc_seq(base_qualities, "r")		
	else:
		print >> fastq2, read2.query_sequence
		print >> fastq2, "+"
		print >> fastq2, base_qualities
	#print read1.tostring(star_bam).strip()	
	#print read2.tostring(star_bam).strip()	

parser = argparse.ArgumentParser()
parser.add_argument('star_bam_file', action='store', help='star alignment')
#parser.add_argument('ensbed', action='store', help='gene annotation')
parser.add_argument('out_fastq1', action='store', help='output fastq1 name')
parser.add_argument('out_fastq2', action='store', help='output fastq2 name')

args = parser.parse_args()

star_bam =  pysam.AlignmentFile(args.star_bam_file, "rb")
read1 = ""
read2 = ""
uniq_quality = 255
#gene_ann = pygeneann.GeneAnnotation(args.ensbed)
with open(args.out_fastq1, 'w') as fastq1, open(args.out_fastq2, 'w') as fastq2:
	for read in star_bam.fetch(until_eof=True):
		
		if read.is_secondary or read.is_supplementary:
			continue
		
		if not read1:
			read1 = read
		elif not read2:
			read2 = read
			if read1.query_name != read2.query_name:
				raise Exception("Singleton read found:" + read1.query_name)
			if read1.mapping_quality != uniq_quality or read2.mapping_quality != uniq_quality:
				#output_paired_fastq(read1, read2, outfq1, outfq2)
				output_paired_fastq(read1, read2, star_bam, fastq1, fastq2)
			elif read1.reference_name != read2.reference_name:
				output_paired_fastq(read1, read2, star_bam, fastq1, fastq2)
			elif abs(read1.template_length) > 300:
				output_paired_fastq(read1, read2, star_bam, fastq1, fastq2)
			elif (read1.reference_start < read2.reference_start and read1.is_reverse) or (read1.reference_start > read2.reference_start and read2.is_reverse) or (read1.is_reverse == read2.is_reverse):
				output_paired_fastq(read1, read2, star_bam, fastq1, fastq2)
			elif (read1.cigarstring != str(read1.query_length) + "M") or (read2.cigarstring != str(read2.query_length) + "M"):
				output_paired_fastq(read1, read2, star_bam, fastq1, fastq2)
				
				"""
				matched_genes1 = gene_ann.map_pos_to_genes(format_chr(read1.reference_name), read1.reference_start)
				matched_genes2 = gene_ann.map_pos_to_genes(format_chr(read2.reference_name), read2.reference_start)
				gene_names1 = set([g.gene_name for g in matched_genes1])
				gene_names2 = set([g.gene_name for g in matched_genes2])
				if gene_names1 and gene_names1 != gene_names2:
					output_paired_fastq(read1, read2, star_bam, fastq1, fastq2)
				elif not gene_names1 and not gene_names2:
					output_paired_fastq(read1, read2, star_bam, fastq1, fastq2)
				"""

			read1 = ""
			read2 = ""
		
