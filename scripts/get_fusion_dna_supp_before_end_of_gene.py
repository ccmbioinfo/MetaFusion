#!/usr/bin/env python

import sys
import os
import commands
import argparse
import pysam

#sys.path.append("/hpf/largeprojects/ccmbio/jiangyue/DIPG_analysis_by_samples/Scripts/pygeneann/pygenefusionann")
import pygeneann

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

parser = argparse.ArgumentParser()

parser.add_argument('cff_file', action='store', help='CFF file (.cff or .cff.reann)')
parser.add_argument('dna_bam_list', action='store', help='A file include paths of DNA sequencing files')
parser.add_argument('gene_bed', action='store', help='Gene annotation bed file')
parser.add_argument('tmp_dir', action='store', help='Temp file directory')

args = parser.parse_args()



cff_file = args.cff_file
dna_bam_list = args.dna_bam_list
out_dir = args.tmp_dir
gene_ann_file = args.gene_bed

#dna_bam =  pysam.AlignmentFile(dna_bam_file, "rb")
bam_dict = load_bam_dict(dna_bam_list)

gene_ann = pygeneann.GeneAnnotation(gene_ann_file)


win_size = 100000
rlen = 100
isize = 500
supp_cluster_num = -99 # No dna file: -1; gene not in annotation: -2; chr not in bam's refrerence: -3; confilicting windonw start and end: -4 
for line in open(cff_file, "r"):
	fusion = pygeneann.CffFusion(line)
	if fusion.sample_name in bam_dict:
		dna_bam_file = bam_dict[fusion.sample_name]
		dna_bam =  pysam.AlignmentFile(dna_bam_file, "rb")
	else:
		#print >> sys.stderr, fusion.sample_name, "has no dna file."
		supp_cluster_num = -1 
		#print "Fusion:", line.strip(), supp_cluster_num
		fusion.dnasupp = supp_cluster_num
		print fusion.tostring()
		continue
	discordant_bam = pysam.Samfile("discordant.bam", "wb", template=dna_bam)

	# build window to extract fusion supporting reads
	interval_t_gene_1 = gene_ann.get_gene_interval(fusion.t_gene1)
	interval_t_gene_2 = gene_ann.get_gene_interval(fusion.t_gene2)
	interval_reann_gene_1 = gene_ann.get_gene_interval(fusion.reann_gene1)
	interval_reann_gene_2 = gene_ann.get_gene_interval(fusion.reann_gene2)
	if interval_t_gene_1.overlap(interval_reann_gene_1):
		interval_1 = interval_t_gene_1.merge(interval_reann_gene_1)
	else:
		interval_1 = interval_reann_gene_1
	if interval_t_gene_2.overlap(interval_reann_gene_2):
		interval_2 = interval_t_gene_2.merge(interval_reann_gene_2)
	else:
		interval_2 = interval_reann_gene_2
		
	# Gene not in annotation
	if interval_1.gene_name == "NA" or interval_2.gene_name == "NA":
		supp_cluster_num = -2
		#print "Fusion:", line.strip(), supp_cluster_num
		fusion.dnasupp = supp_cluster_num
		print fusion.tostring()
		continue

        # FORMAT CHROMOSOMES
	win_chr1 = fusion.chr1[3:]
	win_chr2 = fusion.chr2[3:]
	
	win_strand_reverse1 = False if fusion.strand1 == "+" else True
	win_strand_reverse2 = False if fusion.strand2 == "+" else True
	if fusion.strand1 == "+" and fusion.strand2 == "+":
		win_start1 = fusion.pos1 - rlen
		#max next exon start
		win_end1 = interval_1.end
		#win_end1 = fusion.pos1 + win_size


		win_start2 = fusion.pos2 - rlen
		#win_end2 = fusion.pos2 + win_size
		#max previous exon start
		win_end2 = interval_2.end
		
	elif fusion.strand1 == "+" and fusion.strand2 == "-":
		win_start1 = fusion.pos1 - rlen
		#win_end1 = fusion.pos1 + win_size
		win_end1 = interval_1.end

		#win_start2 = fusion.pos2 - win_size
		win_start2 = interval_2.start
		win_end2 = fusion.pos2 + rlen 
	
	elif fusion.strand1 == "-" and fusion.strand2 == "+":
		#win_start1 = fusion.pos1 - win_size
		win_start1 = interval_1.start
		win_end1 = fusion.pos1 + rlen

		win_start2 = fusion.pos2 - rlen
		#win_end2 = fusion.pos2 + win_size 
		win_end2 = interval_2.end
	elif fusion.strand1 == "-" and fusion.strand2 == "-":
		#win_start1 = fusion.pos1 - win_size
		win_start1 = interval_1.start
		win_end1 = fusion.pos1 + rlen

		#win_start2 = fusion.pos2 - win_size
		win_start2 = interval_2.start
		win_end2 = fusion.pos2 + rlen
	else:
		#print >> sys.stderr, fusion.strand1, fusion.strand2
		#sys.exit(1)
		continue
	
	print >> sys.stderr,"Window:", win_chr1, win_start1, win_end1, win_chr2, win_start2, win_end2
	if dna_bam.get_tid(win_chr1) < 0 or dna_bam.get_tid(win_chr2) < 0: 
		#print "Fusion:", line.strip(), "No_data"
		supp_cluster_num = -3
		#print "Fusion:", line.strip(), supp_cluster_num
		fusion.dnasupp = supp_cluster_num
		print fusion.tostring()
		continue	

	clusters = []
	if win_start1 > win_end1 or win_start2 > win_end2:
		print >> sys.stderr, win_chr1, win_start1, win_end1, win_chr2, win_start2, win_end2
		print >> sys.stderr, line.strip()
		supp_cluster_num = -4
		#print "Fusion:", line.strip(), supp_cluster_num
		fusion.dnasupp = supp_cluster_num
		print fusion.tostring()
		continue
	win_start1 = max(win_start1, 0)
        win_start2 = max(win_start2, 0)
	for read in dna_bam.fetch(win_chr1, win_start1, win_end1):
		if read.is_unmapped or read.mate_is_unmapped:
			continue
		else:
			# if both reads from the same chr, read1.pos should be less than read2.pos, in case window1 and window2 overlap
			read_win_start2 = win_start2
			read_win_end2 = win_end2
			if win_chr1 == win_chr2:
				if read.reference_start < read.next_reference_start:
					read_win_start2 = max(win_start2, read.reference_start)
				else:
					read_win_end2 = min(win_end2, read.reference_start)
			# debug
			# read1 in window1, read2 in window2, strand consistant with fusion strands
			#if read.query_name == "HWI-ST915:90:D04LBACXX:6:1201:12067:170887":
				#print "Found", read.reference_start > win_start1, read.is_reverse == win_strand_reverse1, read.next_reference_name == win_chr2, win_start2 < read.next_reference_start < win_end2, read.mate_is_reverse == win_strand_reverse2 
			#	print win_start2, read.reference_start, read.next_reference_start, win_end2

			if read.reference_start > win_start1 and read.is_reverse == win_strand_reverse1 and read.next_reference_name == win_chr2 and read_win_start2 < read.next_reference_start < read_win_end2 and read.mate_is_reverse == win_strand_reverse2 and (win_chr1 != win_chr2 or abs(read.template_length) > 1000):
				#discordant_bam.write(read)

				# cluster discordant pairs
				#discordant_bam = pysam.Samfile("discordant.bam", "rb")
				chr1 = read.reference_name
				pos1 = read.reference_start
				chr2 = read.next_reference_name
				pos2 = read.next_reference_start
				clustered = False
				for cluster  in clusters:
					if chr1 == cluster.chr1 and chr2 == cluster.chr2 and abs(pos1 - cluster.pos1) < rlen and abs(pos2 - cluster.pos2) < isize:
						cluster.pos1 = (pos1 + cluster.pos1)/2
						cluster.pos2 = (pos2 + cluster.pos2)/2
						cluster.number += 1
						clustered = True
						break
				if not clustered:
					cluster = Cluster(chr1, pos1, chr2, pos2, 1)
					clusters.append(cluster)
				
			
	discordant_bam.close()	
	pysam.sort("discordant.bam", out_dir+"/discordant")
	pysam.index("discordant.bam")

	supp_cluster_num = 0
	for cluster in clusters:
		if cluster.number > 2:
			print >> sys.stderr, "Cluster:", cluster.tostring()
			supp_cluster_num += 1
	#--- output dna supp fusion only ---#
	#if supp_cluster_num:
	#	print "Fusion:", line.strip(), supp_cluster_num

	#--- output all fusions ---#
	fusion.dnasupp = supp_cluster_num
	print fusion.tostring()
	#print "Fusion:", line.strip(), supp_cluster_num
	#print >> sys.stderr, script_dir + "/cluster2.sh " + out_dir + "/discordant.bam header " + out_dir
	#os.system(script_dir + "/cluster2.sh " + out_dir + "/discordant.bam header " + out_dir)

	#a,b = commands.getstatusoutput("wc -l " + out_dir + '/tmp_bps')
	#os.system("cat " + out_dir + "/tmp_bps")
	#tmp = line.split()
	#print "Fusion:", line.strip(), b.split()[0]

	dna_bam.close()
