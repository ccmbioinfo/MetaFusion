#!/usr/bin/env python

import sys
import pysam
import bisect
import time
import sequtils
import copy
import re
class GeneBed():
	def __init__(self, bed_line):
		if bed_line:
			tmp = bed_line.split()
			self.chr = tmp[0]
			self.start = int(tmp[1]) + 1	# bed file use 0-based coordinate
			self.end = int(tmp[2])		# start and end are first and last base of each segment
			self.transcript_id = tmp[3]
			self.type = tmp[4] # utr/intron/cds
			self.idx = int(tmp[5])
			self.strand = tmp[6] # r/f
			self.gene_name = tmp[7]
			self.gene_id = tmp[8]
			# the following two attr are only used for breakpoint annotation, tell whether current breakpoint on/close to boundary, need to reset everytime before used
			#self.is_on_boundary = False
			#self.close_to_boundary = False
			# in gene intervals, this gene includes annotations from different chr, strand or >5Mb distance
			self.is_from_contradictory_gene = False
		else:
			self.chr = "NA"
			self.start = 0	# bed file use 0-based coordinate
			self.end = 0		# start and end are first and last base of each segment
			self.transcript_id = "NA"
			self.type = "NA" # utr/intron/cds
			self.idx = -1
			self.strand = "NA"
			self.gene_name = "NA"
			self.gene_id = "NA"
			self.is_from_contradictory_gene = False

	def overlap(self, genebed2):
		if self.chr == genebed2.chr and min(self.end, genebed2.end) - max(self.start, genebed2.start) > 0:
			return True
		else:
			return False
	def merge(self, genebed2):
		if self.chr != genebed2.chr or self.strand != genebed2.strand:
			print >> sys.stderr, "Warning: intervals are on different chr/strand."
			print >> sys.stderr, self.gene_name, genebed2.gene_name
			return self
		else:
			#new_interval = GeneIntervals("Merged_" + self.gene_name + "_" + interval2.gene_name, self.chr, min(self,start, interval2.start), max(self.end, interval2.end), self.strand, False)
			new_bed = copy.copy(genebed2)
			new_bed.transcript_id = "Merged"
			new_bed.start =  min(self.start, genebed2.start)
			new_bed.end =  max(self.end, genebed2.end)
			
			return new_bed
	def tostring(self):
		attrs = [self.chr, str(self.start), str(self.end), self.transcript_id, self.type, str(self.idx), self.strand, self.gene_name, self.gene_id]
		return "\t".join(attrs)
# fusions are clusted by gene pairs
#Gene_Cluster DIPG106N Normal DIPG EricScript GeneFusion NoDriverGene,GeneFusion >>CTD-2369P2.8_utr3;ICAM1_utr3>>VPS39_intron,utr3 False False False False -1; # old category format

#Gene_Cluster S1PR5 KEAP1 Tumor DIPG Integrate,EricScript ReadThrough True True True True 0 DIPG28T,SJHGG143_A,SJHGG075_A; current format

class CategoryFusions():
	def __init__(self, category_line):
		self.__load_category(category_line)
	def __load_category(self, category_line):
		tmp = category_line.split()
		self.cluster_type = tmp[0]
		self.gene1 = tmp[1]
		self.gene2 = tmp[2]
		self.sample_type = tmp[3] #tmp[2].split(",")
		self.disease = tmp[4].split(",")
		self.tools = tmp[5].split(",")
		self.inferred_fusion_type = tmp[6]
		self.gene1_on_bnd = tmp[7]
		self.gene1_close_to_bnd = tmp[8]
		self.gene2_on_bnd = tmp[9]
		self.gene2_close_to_bnd = tmp[10]
		self.dna_supp = tmp[11] # dna support, 0: has data, no dna pairs; >0:  number of dna pair clusters; -1: no data; -2: gene not in annotation; -3: chr not in bam; -4: confilicting windown start and end
		self.samples = tmp[12].split(",")
		self.line = category_line.strip()
	def out(self):
		print self.line
		
			
class CategoryFusionStats():
	category_list = []
	def __init__(self, category_file):
		self.__load_category_file(category_file)
	def __load_category_file(self, category_file):
		for line in open(category_file, "r"):
			if line.startswith("#"):
				continue
			category = CategoryFusions(line)
			self.category_list.append(category)
	def filter_recurrent(self, fusion_list, threshold):
		return filter(lambda x:len(x.samples)>=threshold, fusion_list)
	def filter_disease(self, fusion_list, disease):
		return filter(lambda x:disease in x.disease, fusion_list)

	def filter_sample_number(self, fusion_list, threshold, sample_prefix):
		return filter(lambda x:len(filter(lambda y:y.startswith(sample_prefix), x.samples)) >= threshold, fusion_list)

	def filter_sample_type (self, fusion_list, sample_type): 	
		return filter(lambda x:x.sample_type==sample_type, fusion_list)
	def filter_tools_name (self, fusion_list, tool_name): 	
		return filter(lambda x:tool_name in x.tools, fusion_list)
	def filter_tools_num (self, fusion_list, tools_num): 	
		return filter(lambda x:len(x.tools)>=tools_num, fusion_list)
	def filter_inferred_type (self, fusion_list, inferred_type): 	
		return filter(lambda x:x.inferred_fusion_type==inferred_type, fusion_list)
	def filter_on_bnd (self, fusion_list): 	
		return filter(lambda x:x.gene1_on_bnd and x.gene2_on_bnd, fusion_list)
	def filter_close_to_bnd (self, fusion_list): 	
		return filter(lambda x:(x.gene1_close_to_bnd == "True" and x.gene2_close_to_bnd == "True"), fusion_list)
	def filter_dna_supp (self, fusion_list): 	
		return filter(lambda x:int(x.dna_supp)>0, fusion_list)
		
	


			
class CffFusionStats():
	__fusion_dict = {}
	__fusion_samples_dict = {}
	def __init__(self, cff_file):
		pass
		#self.__load_fusions(cff_file)
	def __load_fusions(self, cff_file):
		for line in open(cff_file, "r"):
			fusion = CffFusion(line)
			key = fusion.t_gene1 + "_" + fusion.t_gene2
			rkey = fusion.t_gene2 + "_" + fusion.t_gene1
			if rkey not in self.__fusion_dict:
				self.__fusion_dict.setdefault(key, []).append(fusion)
				self.__fusion_samples_dict.setdefault(key, []).append(fusion.sample_name)
			else:
				self.__fusion_dict.setdefault(rkey,[]).append(fusion)
				self.__fusion_samples_dict.setdefault(rkey, []).append(fusion.sample_name)
	
	# compare two fusions based on the up/downstream gene sets, if overlap on both sets consider them as same fusion
	def is_same_gene_pair_fusion(self, fusion1, fusion2):
		id = 1 
		# compare gene sets on fw strand 
		up_g1, down_g1 = set(fusion1.get_reannotated_genes(id))
		up_g2, down_g2 = set(fusion2.get_reannotated_genes(id))

		if up_g1 & up_g2 and down_g1 & down_g2:
			return True

		id = 2 
		# compare gene sets on bw strand 
		up_g1, down_g1 = set(fusion1.get_reannotated_genes(id))
		up_g2, down_g2 = set(fusion2.get_reannotated_genes(id))

		if up_g1 & up_g2 and down_g1 & down_g2:
			return True

		return False
	# compare gene fusions based on breakpoints
	def generate_common_fusion_stats(self, cff_file):
		fusion_bp_dict = {}
		fusion_pos1_idx_dict = {}
		# build index for fusions based on fusion.pos1
		for line in  open(cff_file, "r"):
			if line.startswith("#"):
				continue
			fusion = CffFusion(line)
			fusion_bp_dict.setdefault(fusion.chr1, {}).setdefault(fusion.strand1, []).append(fusion)
		for chr in fusion_bp_dict:
			for strand in fusion_bp_dict[chr]:
				fusion_bp_dict[chr][strand].sort(key = lambda x:x.pos1)	
				fusion_pos1_idx_dict.setdefault(chr, {}).setdefault(strand, [f.pos1 for f in fusion_bp_dict[chr][strand]]) 
		# find common fusions
	# output fusions in a fusion list as clustered fusions
	def output_clustered_fusions(self, fusion_list, cluster_type):	
			gene1_list = [f.reann_gene1 for f in fusion_list]
			gene2_list = [f.reann_gene2 for f in fusion_list]

			sample_list = [f.sample_name for f in fusion_list]	
			disease_list = [f.disease for f in fusion_list]	
			tool_list = [f.tool for f in fusion_list]	
			sample_type_list = [f.sample_type for f in fusion_list]	
			gene1_on_bndry = "True" in [f.gene1_on_bndry for f in fusion_list] 
			gene1_close_to_bndry = "True" in [f.gene1_close_to_bndry for f in fusion_list] 
			gene2_on_bndry = "True" in [f.gene2_on_bndry for f in fusion_list] 
			gene2_close_to_bndry = "True" in [f.gene2_close_to_bndry for f in fusion_list] 
			
			dna_supp_cluster_num = max([int(f.dnasupp) for f in fusion_list])

			category_list = [f.category for f in fusion_list]	

			print cluster_type, ",".join(list(set(gene1_list))), ",".join(list(set(gene2_list))), ",".join(list(set(sample_type_list))), ",".join(list(set(disease_list))), ",".join(list(set(tool_list))), ",".join(list(set(category_list))), gene1_on_bndry, gene1_close_to_bndry, gene2_on_bndry, gene2_close_to_bndry, dna_supp_cluster_num, ",".join(list(set(sample_list)))

	def output_clustered_fusions_version1(self, fusion_list, cluster_type):	
			sample_list = [f.sample_name for f in fusion_list]	
			disease_list = [f.disease for f in fusion_list]	
			tool_list = [f.tool for f in fusion_list]	
			sample_type_list = [f.sample_type for f in fusion_list]	
			gene1_on_bndry = "True" in [f.gene1_on_bndry for f in fusion_list] 
			gene1_close_to_bndry = "True" in [f.gene1_close_to_bndry for f in fusion_list] 
			gene2_on_bndry = "True" in [f.gene2_on_bndry for f in fusion_list] 
			gene2_close_to_bndry = "True" in [f.gene2_close_to_bndry for f in fusion_list] 
			
			dna_supp_cluster_num = max([int(f.dnasupp) for f in fusion_list])



			category_list = [f.reann_category1 for f in fusion_list]	
			category_list += [f.reann_category2 for f in fusion_list]	
			gene_order_list = [f.reann_gene_order1 for f in fusion_list]
			gene_order_list += [f.reann_gene_order2 for f in fusion_list]
			
			inferred_category = "None"
			# infer possible gene fusion category
			for category in ["SameGene", "ReadThrough", "GeneFusion", "TruncatedCoding", "TruncatedNoncoding", "NoDriverGene"]:
				if category in category_list:
					inferred_category = category
					break

			print cluster_type, ",".join(list(set(sample_list))), ",".join(list(set(sample_type_list))), ",".join(list(set(disease_list))), ",".join(list(set(tool_list))), inferred_category, ",".join(list(set(category_list))), "|".join(list(set(gene_order_list))), gene1_on_bndry, gene1_close_to_bndry, gene2_on_bndry, gene2_close_to_bndry, dna_supp_cluster_num
	# cluster fusions of "NoDriverGene" and "Truncated" type on their breakpoints
	def generate_common_fusion_stats_by_breakpoints(self, fusion_list):
		diff = 100000
		# save clustered fusion id, skip a fusion when it is already clustered
		clustered_id = {}
		print >> sys.stderr, "bp fusion list:", len(fusion_list)
		for i in range(len(fusion_list)):
			if i in clustered_id:
				continue

			#if i % 100 == 0: print >> sys.stderr, i

			fusion1 = fusion_list[i]
			small_bp1, big_bp1 = fusion1.get_ordered_breakpoints()
			clustered_id.setdefault(i, i)
			fusion_cluster_list = []
			fusion_cluster_list.append(fusion1)
			for j in range(len(fusion_list)):
				if j in clustered_id:
					continue
				if i == j:
					continue

				fusion2 = fusion_list[j]
				small_bp2, big_bp2 = fusion2.get_ordered_breakpoints()
				if cmp_fusion_breakpoints(small_bp1, small_bp2, diff) and cmp_fusion_breakpoints(big_bp1, big_bp2, diff):
					clustered_id.setdefault(j, j)
					fusion_cluster_list.append(fusion2)

			self.output_clustered_fusions(fusion_cluster_list, "BP_Cluster")
				

	# compare gene fusions based on re-annotated gene names, hard to work for truncated fusions, used generate_common_fusion_stats_by_breakpoints for NoDriverGenes and Trucated type fusions
	def generate_common_fusion_stats_by_genes(self, cff_file):
		fusion_dict = {}
		fusion_list_for_bp_cmp = []
		common_key_dict = {}
		# cluster fusions by gene pairs, save in fusion_dict
		for line in open(cff_file, "r"):
			if line.startswith("#"):
				continue
			fusion = CffFusion(line)
			# send fusion to breakpoint cluster later
			if fusion.reann_gene1 == "NA" or fusion.reann_gene2 == "NA":
				fusion_list_for_bp_cmp.append(fusion)
			else:
				key = (fusion.reann_gene1, fusion.reann_gene2)	
				fusion_dict.setdefault(key, []).append(fusion)
		# output clustered fusions
		for key in fusion_dict:
			fusion_list = fusion_dict[key]

			self.output_clustered_fusions(fusion_list, "Gene_Cluster")
		self.generate_common_fusion_stats_by_breakpoints(fusion_list_for_bp_cmp)		
		
	def generate_common_fusion_stats_by_genes_version1(self, cff_file):
		fusion_dict = {}
		fusion_list_for_bp_cmp = []
		common_key_dict = {}
		for line in open(cff_file, "r"):
			if line.startswith("#"):
				continue
			fusion = CffFusion(line)
			keys = set()
			found = False # found both up/down stream genes
			id = 1
			if fusion.reann_category1 not in ["SameGene", "NoDriverGene"]:
				up_g, down_g = fusion.get_reannotated_genes(id)
				for g1 in up_g:
					for g2 in down_g:
						key = g1 + "_" + g2 
						keys.add(key)
						fusion_dict.setdefault(key, []).append(fusion)
						found = True
			id = 2
			if fusion.reann_category2 not in ["SameGene", "NoDriverGene"]:
				up_g, down_g = fusion.get_reannotated_genes(id)
				for g1 in up_g:
					for g2 in down_g:
						key = g1 + "_" + g2 
						keys.add(key)
						fusion_dict.setdefault(key, []).append(fusion)
						found = True
			# can not find a meaningful gene pair for this fuison, try breakpoints comparison later
			if not found:
				fusion_list_for_bp_cmp.append(fusion)
			# for breakpoints that can be mapped to multiple genes pairs, save these pairs in a dict, merge their fusion lists when output
			if len(keys) > 1:
				for key in keys:
					if key in common_key_dict:
						common_key_dict[key] |= keys
					else:
						common_key_dict[key] = keys
		removed_keys = {} 
		for key in fusion_dict:
			if key in removed_keys:
				continue
			fusion_list = fusion_dict[key]

			if key in common_key_dict:
				key_set = common_key_dict[key]
				while True:
					l = len(key_set)
					key_set2 = key_set.copy()
					for key2 in key_set2:
						if key2 != key:
							key_set |= common_key_dict[key2]	
					if l == len(key_set):
						break
				#print key, key_set
				key_set.remove(key)
				for key2 in key_set:
					fusion_list += fusion_dict[key2]
					removed_keys.setdefault(key2, key2)	
			self.output_clustered_fusions(fusion_list, "Gene_Cluster")
			'''
			sample_list = [f.sample_name for f in fusion_list]	
			disease_list = [f.disease for f in fusion_list]	
			tool_list = [f.tool for f in fusion_list]	
			sample_type_list = [f.sample_type for f in fusion_list]	

			category_list = [f.reann_category1 for f in fusion_list]	
			category_list += [f.reann_category2 for f in fusion_list]	
			gene_order_list = [f.reann_gene_order1 for f in fusion_list]
			gene_order_list += [f.reann_gene_order2 for f in fusion_list]
			print key, ",".join(list(set(sample_list))), ",".join(list(set(sample_type_list))), ",".join(list(set(disease_list))), ",".join(list(set(tool_list))), ",".join(list(set(category_list))), ",".join(list(set(gene_order_list)))
			'''
		# send "NoDriverGene" and "Truncated" type fusions for breakpoint cluster
		self.generate_common_fusion_stats_by_breakpoints(fusion_list_for_bp_cmp)		
	def get_gene_order_stats(self):
		for key in self.__fusion_dict:
			n_sg = 0 # same gene
			n_rt = 0 # read through
			n_gf = 0 # gene fusion
			n_tc = 0 # truncated coding
			n_tn = 0 # truncated noncoding
			n_ns = 0 # nonsense
			gene_order = []
			sample_type = []
			for sample in set(self.__fusion_samples_dict[key]):
				if sample.endswith("N"):
					sample_type.append("Normal")
				elif sample.endswith("T"):
					sample_type.append("Tumor")
				else:
					sample_type.append("Unknown")

			type = "Unknown"
			for fusion in self.__fusion_dict[key]:
				#tmp = ";".join(fusion.otherann)
				tmp = []
				for attr in fusion.zone4_attrs:
					tmp.append(fusion.__dict__[attr])
				gene_order.append(" ".join(tmp))
				if "SameGene" in tmp:
					type = "SameGene"
					n_sg += 1
				elif "ReadThrough" in tmp:
					type = "ReadThrough" if type == "Unknown" else type
					n_rt += 1
				elif "GeneFusion" in tmp:
					type = "GeneFusion" if type == "Unknown" else type
					n_gf += 1
				elif "TruncatedCoding" in tmp:
					type = "TruncatedCoding" if type == "Unknown" else type
					n_tc +=1
				elif "TruncatedNoncoding" in tmp:
					type = "TruncatedNoncoding" if type == "Unknown" else type
					n_tn += 1
				elif "NoDriverGene" in tmp:
					type = "NoDriverGene" if type == "Unknown" else type
					n_ns += 1
				else:
					print >> sys.stderr, "Fusions without category:", fusion.tostring()
			if type != "Unknown":
				print "Fusion", fusion.t_gene1, fusion.t_gene2, ",".join(list(set(self.__fusion_samples_dict[key]))), type, ",".join(sorted(list(set(sample_type)))),
				print "|".join(list(set(gene_order)))
				print "\tSameGene:", n_sg
				print "\tReadThrough:", n_rt
				print "\tGeneFusion:", n_gf
				print "\tTruncatedCoding:", n_tc
				print "\tTruncatedNoncoding:", n_tn
				print "\tNonSense:", n_ns


class CffFusion():
	# chr1 pos1 strand1 chr2 pos2 strand2 library sample_name sample_type disease tool split_cnt span_cnt tool_gene1 tool_gene_area1 tool_gene2 tool_gene_area2 fwd_gene_order bwd_gene_order
	def __init__(self, cff_line):
		tmp = cff_line.split()
		# Breadkpoint Zone
		self.chr1 = tmp[0]
		self.pos1 = int(tmp[1])
		self.strand1 = tmp[2]
		self.chr2 = tmp[3]
		self.pos2 = int(tmp[4])
		self.strand2 = tmp[5]
		# Sample info Zone
		self.library = tmp[6] # DNA/RNA
		self.sample_name = tmp[7]
		self.sample_type = tmp[8] # Tumor/Normal
		self.disease = tmp[9]
		# Software Zone
		self.tool = tmp[10]
		self.split_cnt = int(tmp[11])
		self.span_cnt = int(tmp[12]) if tmp[12] != "NA" else tmp[12]
		self.t_gene1 = tmp[13] # gene reported by tool
		self.t_area1 = tmp[14] # exon/utr/intron
		self.t_gene2 = tmp[15]
		self.t_area2 = tmp[16]
		# Re-annotation Zone
		# ReadThrough     DTX2    cds     DTX2P1-UPK3BP1-PMS2P11  utr3    True    TrueTrue     True    5.5     1       474827  1       F00000001       CCTCCCGCAGGGCCCTGAGCACCCCAATCCCGGAAAGCCGTTCACTGCCAGAGGGTTTCCCCGCCAGTGCTACCTTCCAGACAACGCCCAGGGCCGCAAG    CCTCCAGGGGCTTCCAGAACCCGGAGACACTGGCTGACATTCCGGCCTCCCCACAGCTGCTGACCGATGGCCACTACATGACGCTGCCCGTGTCTCCGGA
		#if len(tmp) == 33:
		if len(tmp) >= 33 :
			self.category = tmp[17]
			self.reann_gene1 = tmp[18]
			self.reann_type1 = tmp[19]
			self.reann_gene2 = tmp[20]
			self.reann_type2 = tmp[21]
			self.gene1_on_bndry = tmp[22] # gene1 breakpoint on boundary
			self.gene1_close_to_bndry = tmp[23] # gene1 breakpoint within 10bp of a boundary
			self.gene2_on_bndry = tmp[24] # gene2 on boundary
			self.gene2_close_to_bndry = tmp[25] # gene2 close to boundary
			self.score = tmp[26]	# score for bpann1 + bpann2
			self.coding_id_distance = tmp[27]	# difference between two coding fusion genes, if not coding the value is -1
			self.gene_interval_distance = tmp[28] # distance between two fusion gene intervals
			self.dnasupp = tmp[29]
			self.fusion_id = tmp[30]
			self.seq1 = tmp[31]
			self.seq2 = tmp[32]
			self.is_inframe = False
			self.splice_site1 = "NA"
			self.splice_site2 = "NA"
		else:
			self.category = "NA"	# category
			self.reann_gene1 = "NA"
			self.reann_type1 = "NA"
			self.reann_gene2 = "NA"
			self.reann_type2 = "NA"
			self.gene1_on_bndry = "NA" # gene1 breakpoint on boundary
			self.gene1_close_to_bndry = "NA" # gene1 breakpoint within 10bp of a boundary
			self.gene2_on_bndry = "NA" # gene2 on boundary
			self.gene2_close_to_bndry = "NA" # gene2 close to boundary
			self.score = 0	# score for bpann1 + bpann2
			self.coding_id_distance = -1	# difference between two coding fusion genes, if not coding the value is -1
			self.gene_interval_distance = 0 # distance between two fusion gene intervals
			self.fusion_id = "NA"
			self.seq1 = "NA"
			self.seq2 = "NA"
			self.is_inframe = False
			self.splice_site1 = "NA"
			self.splice_site2 = "NA"
			if len(tmp) == 30:
				self.dnasupp = tmp[29]
			else:
				self.dnasupp = -9 # not available
		self.bpann1 = GeneBed("")	# breakpoint1 annotation (GeneBed)
		self.bpann2 = GeneBed("")	# breakpoint2 annotation (GeneBed)	
		self.score1 = 0	# score for bpann1 + bpann2
		self.score2 = 0	# score for bpann1 + bpann2
		
			
		'''
		if len(tmp) >= 25:
			self.reann_gene_order1 = tmp[17] # re-annotated gene order on fw strand e.g. ZNF248_utr3,cds>>RP11-258F22.1_utr3
			self.reann_gene_type1 = tmp[18] # re-annotated gene type e.g. CodingGene>>NoncodingGene
			self.reann_gene_index1 = tmp[19] # gene index for read-through inference e.g. 7409_r>>9974_r, only valid for coding gene
			self.reann_category1 = tmp[20] # infered fusion type, ReadThrough, GeneFusion, TruncatedCoding, TruncatedNoncoding, Nosenes, SameGene

			self.reann_gene_order2 = tmp[21] # backward re-annotation
			self.reann_gene_type2 = tmp[22] 
			self.reann_gene_index2 = tmp[23] 
			self.reann_category2 = tmp[24]

		else:
			self.reann_gene_order1 = "NA" # re-annotated gene order on fw strand e.g. ZNF248_utr3,cds>>RP11-258F22.1_utr3
			self.reann_gene_type1 = "NA" # re-annotated gene type e.g. CodingGene>>NoncodingGene
			self.reann_gene_index1 = "NA" # gene index for read-through inference e.g. 7409_r>>9974_r, only valid for coding gene
			self.reann_category1 = "NA" # infered fusion type, ReadThrough, GeneFusion, TruncatedCoding, TruncatedNoncoding, Nosenes, SameGene

			self.reann_gene_order2 = "NA" # backward re-annotation
			self.reann_gene_type2 = "NA"
			self.reann_gene_index2 = "NA"
			self.reann_category2 = "NA"
		
		# breakpoints on boundary
		if len(tmp) >= 29 and False:
			self.gene1_on_bndry = tmp[25] # gene1 breakpoint on boundary
			self.gene1_close_to_bndry = tmp[26] # gene1 breakpoint within 10bp of a boundary
			self.gene2_on_bndry = tmp[27] # gene2 on boundary
			self.gene2_close_to_bndry = tmp[28] # gene2 close to boundary
		else:
			self.gene1_on_bndry = "NA" # gene1 breakpoint on boundary
			self.gene1_close_to_bndry = "NA" # gene1 breakpoint within 10bp of a boundary
			self.gene2_on_bndry = "NA" # gene2 on boundary
			self.gene2_close_to_bndry = "NA" # gene2 close to boundary
		# column 30 is dna supporting cluster number
		'''

		self.boundary_info = ""		
		# same all attrs in a list, for printing	
		self.zone1_attrs = ["chr1", "pos1", "strand1", "chr2", "pos2", "strand2"]
		self.zone2_attrs = ["library", "sample_name", "sample_type", "disease"]
		self.zone3_attrs = ["tool", "split_cnt", "span_cnt", "t_gene1", "t_area1", "t_gene2", "t_area2"]
		self.zone4_attrs = ["reann_gene_order1", "reann_gene_type1", "reann_gene_index1", "reann_category1", "reann_gene_order2", "reann_gene_type2", "reann_gene_index2", "reann_category2"]
		# format chr
		if not self.chr1.startswith("chr"):
			self.chr1 = "chr" + self.chr1
		if not self.chr2.startswith("chr"):
			self.chr2 = "chr" + self.chr2
		# check fields
		if not self.check_cff():
			sys.exit(1)
	def get_gene_names_from_gene_order(self):
		# gene name list
		g1 = []
		g2 = []
		g3 = []
		g4 = []
		if self.reann_gene_order1 != "NA":
			tmp = self.reann_gene_order1.split(">>")
			tmp2 = tmp[0].split(";")
			for t2 in tmp2:
				tmp3 = t2.split("_")
				g1.append(tmp3[0])
			tmp2 = tmp[1].split(";")
			for t2 in tmp2:
				tmp3 = t2.split("_")
				g2.append(tmp3[0])
			
			
		if self.reann_gene_order2 != "NA":
			tmp = self.reann_gene_order2.split(">>")
			tmp2 = tmp[0].split(";")
			for t2 in tmp2:
				tmp3 = t2.split("_")
				g3.append(tmp3[0])
			tmp2 = tmp[1].split(";")
			for t2 in tmp2:
				tmp3 = t2.split("_")
				g4.append(tmp3[0])
		
		return g1, g2, g3, g4		

	def __test_prefered_gene(self, genes, gene_name_list):
		prefered_gene = ""
		for gene in genes:
			if gene.gene_name in gene_name_list:
				if not prefered_gene:
					prefered_gene = gene
				else:
					if gene.type == "cds":
						prefered_gene = gene
						break
					elif gene.type == "intron" and prefered_gene != "cds":
						prefered_gene = gene
						
		return prefered_gene		
			
	# among all the re-annotated genes, return a prefered one which should be: 1. breakpoint in/on exon. 2. if not, has a exon closest to breakpoint. 3. the longest isoform.
	# Not finished, needs to consider
	def get_perfered_ann(self, gene_ann):
		print >> sys.stderr,  "Not finished."
		return
		genes1 = gene_ann.map_pos_to_genes(self.chr1, self.pos1)
		genes2 = gene_ann.map_pos_to_genes(self.chr2, self.pos2)
		# get gene name list from re-annotated gene orders
		g1, g2, g3, g4 = self.get_gene_names_from_gene_order()
		
		prefered_gene1 = self.__test_prefered_gene(genes1, g1)
		prefered_gene2 = self.__test_prefered_gene(genes2, g2)
		# try the opposite, because breakpoints order are not always consistant with gene orders
		if not (prefered_gene1 and prefered_gene2):
			prefered_gene1 = self.__test_prefered_gene(genes1, g2)
			prefered_gene2 = self.__test_prefered_gene(genes2, g1)
			
		prefered_gene3 = self.__test_prefered_gene(genes1, g3)
		prefered_gene4 = self.__test_prefered_gene(genes2, g4)
			
		if not (prefered_gene3 and prefered_gene4):
			prefered_gene3 = self.__test_prefered_gene(genes1, g4)
			prefered_gene4 = self.__test_prefered_gene(genes2, g3)
		
		print self.tostring()
		if prefered_gene1 and prefered_gene2:
			print "\t" + prefered_gene1.tostring()
			print "\t" + prefered_gene2.tostring()
		elif prefered_gene3 and prefered_gene4:		
			print "\t" + prefered_gene3.tostring()
			print "\t" + prefered_gene4.tostring()
		else:
			print "No prefered genes."
		
	# compare fusion breakpoints and return in an order of smaller to  bigger
	def get_ordered_breakpoints(self):
		if self.chr1 < self.chr1:
			small_bp = (self.chr1, self.pos1, self.strand1)
			big_bp = (self.chr2, self.pos2, self.strand2)
		elif self.chr1 == self.chr1 and self.pos1 < self.pos2:
			small_bp = (self.chr1, self.pos1, self.strand1)
			big_bp = (self.chr2, self.pos2, self.strand2)
		else:
			small_bp = (self.chr2, self.pos2, self.strand2)
			big_bp = (self.chr1, self.pos1, self.strand1)
		return small_bp, big_bp	
				
	# after reannotation, get a list of upstream genes and a list of downstream genes, can be used to compare fusions from different tools on gene leverl
	def get_reannotated_genes(self, id):
		up_genes = []
		down_genes = []
		if id == 1:
			gene_order = self.reann_gene_order1
		elif id == 2:
			gene_order = self.reann_gene_order2
		else:
			print >> sys.stderr, "Strand has to be 1 or 2", id, "provided."
			sys.exit(1)
			
		#LINC00875_utr5,intron>>NBPF9_intron;LOC100288142_intron;NBPF8_intron    NoncodingGene>>CodingGene,CodingGene,CodingGene >>566_f,525_r,565_f     TruncatedNoncoding
		if gene_order != "NA":
			tmp = gene_order.split(">>")
			for g in tmp[0].split(";"):
				gname = g.split("_")[0]
				if g: up_genes.append(gname)
			for g in tmp[1].split(";"):
				gname = g.split("_")[0]
				if g: down_genes.append(gname)
		'''
		if not up_genes:
			up_genes.append("InterGenic")		
		if not down_genes:
			down_genes.append("InterGenic")		
		'''
		return set(up_genes), set(down_genes)
	def check_cff(self):
		if self.library not in ["NA", "DNA", "RNA"]:
			print >> sys.stderr, "Unknown library type:", self.library
			return False
		if self.sample_type not in ["NA", "Tumor", "Normal"]:
			print >> sys.stderr, "Unknown sample type:", self.sample_type
			return False
		return True

	def tostring(self):
		value = []
		#for attr in self.zone1_attrs + self.zone2_attrs + self.zone3_attrs + self.zone4_attrs:
		for attr in self.zone1_attrs + self.zone2_attrs + self.zone3_attrs:
			if not attr in self.__dict__:
				print >> sys.stderr, "Attribute name error:", attr
				sys.exit(1)
			else:
				value.append(self.__dict__[attr])
		#return "\t".join(map(lambda x:str(x), value)) + "\t" + self.boundary_info
		self.boundary_info = "\t".join(map(str, [self.gene1_on_bndry, self.gene1_close_to_bndry, self.gene2_on_bndry, self.gene2_close_to_bndry]))
		
		return "\t".join(map(lambda x:str(x), value)) + "\t" + self.category + "\t" + self.reann_gene1 + "\t" + self.reann_type1 + "\t" + self.reann_gene2 + "\t" + self.reann_type2 + "\t" + self.boundary_info + "\t" + str(self.score) + "\t" + str(self.coding_id_distance) + "\t" + str(self.gene_interval_distance) + "\t" + str(self.dnasupp) + "\t" + self.fusion_id + "\t" + self.seq1 + "\t" + self.seq2 + "\t" + str(self.is_inframe) + "\t" + self.splice_site1 + "\t" + self.splice_site2
	
	def __check_boundary(self, bpann, order): # bpann is GeneBed object, order = head/tail
		# set on boundary info, used boundaries according to head/tail gene and their strand
		is_on_boundary = False
		close_to_boundary = False
		t = 5
	
		if bpann.start <= self.pos1 <= bpann.end:
			pos = self.pos1
		elif bpann.start <= self.pos2 <= bpann.end:
			pos = self.pos2
		else:
			#print >> sys.stderr, "Map pos to gene error"
			#sys.exit(1)
			return is_on_boundary, close_to_boundary
			
		
		if order == "head":
			if bpann.strand == "f": # + strand head gene
				if bpann.end == pos:
					is_on_boundary = True
					close_to_boundary = True
				elif abs(bpann.end-pos) < t:
					close_to_boundary = True
			elif bpann.strand == "r":
				if bpann.start == pos:
					is_on_boundary = True
					close_to_boundary = True
				elif abs(bpann.start-pos) < t:
					close_to_boundary = True
			else:
				print >> "Unknkown strand:", bpann.strand
				sys.exit(1)
		else:
			if bpann.strand == "f": # + strand head gene
				if bpann.start == pos:
					is_on_boundary = True
					close_to_boundary = True
				elif abs(bpann.start-pos) < t:
					close_to_boundary = True
			elif bpann.strand == "r":
				if bpann.end == pos:
					is_on_boundary = True
					close_to_boundary = True
				elif abs(bpann.end-pos) < t:
					close_to_boundary = True
			else:
				print >> "Unknkown strand:", bpann.strand
				sys.exit(1)
		#print bpann.start, bpann.end, pos
		return is_on_boundary, close_to_boundary
	# assign a score for every potential fusion gene pair;order is head/tail gene		
	def __cal_score(self, bpann, order):
		# Scores to choose best annotation
		score_exon_bnd = 3
		score_exon_bnd_close = 2.9
		score_utr_bnd = 2.5
		score_utr_bnd_close = 2.4
		score_in_exon = 2
		score_in_utr = 1.5
		score_in_intron = 1
		score_intergenic = 0.5
		if bpann.gene_name == "NA":
			score = score_intergenic
			return score
		is_on_boundary, close_to_boundary = self.__check_boundary(bpann, order)
		score = 0
		if bpann.type == "cds":
			if is_on_boundary:
				score = score_exon_bnd
			elif close_to_boundary:
				score = score_exon_bnd_close
			else:
				score = score_in_exon
		elif bpann.type == "utr3" or bpann.type == "utr5":
			if is_on_boundary:
				score = score_utr_bnd
			elif close_to_boundary:
				score = score_utr_bnd_close
			else:
				score = score_in_utr
		elif bpann.type == "intron":
			score = score_in_intron
		else:
			print >> sys.stderr, "Unkknown type:", bpann.type
			sys.exit(1)
		#print score, is_on_boundary, close_to_boundary	
		return score, is_on_boundary, close_to_boundary	
						
	# according to fusion strand (defuse style, strands are supporting pairs') return all possible gene fusions;depending on gene1 and gene2 (a,b,c,d), may have to switch pos order
	def __check_gene_pairs(self, genes1, genes2, gene_ann, switch_pos):
		gene_order = []
		type1 = []
		type2 = []
		id1 = []
		id2 = []
		category = ""
		t = 5


				
		common_genes =  set(genes1.keys()) & set(genes2.keys())

		# type of genes, coding gene ids
		for gene_name in set(genes1.keys()):
			if gene_ann.is_coding(gene_name):
				type1.append("CodingGene")
				id1.append(str(gene_ann.get_coding_gene_idx(gene_name)))
			else:
				type1.append("NoncodingGene")
		for gene_name in set(genes2.keys()):
			if gene_ann.is_coding(gene_name):
				type2.append("CodingGene")
				id2.append(str(gene_ann.get_coding_gene_idx(gene_name)))
			else:
				type2.append("NoncodingGene")
		'''
		# check cds/utr/intron
		#list1 = []
		type1 = ""
		for gene_name in genes1:
			tmp = set([g.type for g in genes1[gene_name]])
			gene1, type1 = self.__check_type(genes1)
			if gene_ann.is_contradictory(gene_name):
				gene_name += "(Cont)"
			#list1.append(gene_name + "_" + ",".join(list(tmp)))
		#list2 = []
		type2 = ""
		for gene_name in genes2:
			tmp = set([g.type for g in genes2[gene_name]])
			gene2, type2 = self.__check_type(genes2)
			if gene_ann.is_contradictory(gene_name):
				gene_name += "(Cont)"
			#list2.append(gene_name + "_" + ",".join(list(tmp)))
		'''
		# No driver gene
		if not genes1:
			category = "NoDriverGene"
		# map to same gene
		elif common_genes:
			category = "SameGene"
		else:
			# category fusions into: read through, gene fusion, truncated coding, truncated noncoding, nonsense
			if "CodingGene" in type1 and "CodingGene" in type2:
				for id in id1:
					tmp = id.split("_")
					idx = int(tmp[0])
					strand = tmp[1]
					#ReadThrough: gene1 and gene2 are adjacent genes (id1 - id2 = 1) or overlapping genes (id1 = id2) but breakpoints cannot map to same gene
					if (strand == "f" and  str(idx+1) + "_f" in id2 ) or (strand == "r" and  str(idx-1) + "_r" in id2) or (id in id2):
						category = "ReadThrough"
					else:
						category = "GeneFusion"
			elif "CodingGene" in type1:
				category = "TruncatedCoding"
			elif "NoncodingGene" in type1:
				category = "TruncatedNoncoding"
			elif not type1:
				category = "NonSense"
			else:
				print >> sys.stderr, "Warning: Unknown category."
				print >> sys.stderr, type1, type2

		# Calsulate score for every gene pair, choose the best
		
		# get best score for current gene combination	
		max_t1 = 0, "NA", "NA", GeneBed("")
		max_t2 = 0, "NA", "NA", GeneBed("")
		for gname1 in genes1:
			for bpann1 in genes1[gname1]:
				score1, is_on_boundary1, close_to_boundary1 = self.__cal_score(bpann1, "head")	
				if score1 > max_t1[0]:
					max_t1 = score1, is_on_boundary1, close_to_boundary1, bpann1
					'''
					self.score1 = score1
					self.gene1_on_bndry = is_on_boundary1
					self.gene1_close_to_bndry = close_to_boundary1
					self.bpann1 = bpann1
					'''
		for gname2 in genes2:
			for bpann2 in genes2[gname2]:
				score2, is_on_boundary2, close_to_boundary2 = self.__cal_score(bpann2, "tail")	
				if score2 > max_t2[0]:
					max_t2 = score2, is_on_boundary2, close_to_boundary2, bpann2
				'''	
				if score2 > self.score2:
					self.score2 = score2
					self.gene2_on_bndry = is_on_boundary2
					self.gene2_close_to_bndry = close_to_boundary2
					self.bpann2 = bpann2
				'''
		if max_t1[0] + max_t2[0] > self.score1 + self.score2:
			self.score1, self.gene1_on_bndry, self.gene1_close_to_bndry, self.bpann1 = max_t1
			self.score2, self.gene2_on_bndry, self.gene2_close_to_bndry, self.bpann2 = max_t2
			self.score = self.score1 + self.score2
			self.category = category
			self.reann_gene1 = self.bpann1.gene_name
			self.reann_gene2 = self.bpann2.gene_name
			self.reann_type1 = self.bpann1.type
			self.reann_type2 = self.bpann2.type

			if id1 and id2:
				if len(set(id1)) != 1 or len(set(id2))!= 1:
					print >> sys.stderr, "coding id err."
					print >> sys.stderr, id1, id2
					sys.exit(1)
				idx1 = int(id1[0].split("_")[0])
				idx2 = int(id2[0].split("_")[0])
				self.coding_id_distance = abs(idx1 - idx2)
			gene_interval1 = gene_ann.get_gene_interval(self.bpann1.gene_name)
			gene_interval2 = gene_ann.get_gene_interval(self.bpann2.gene_name)
			if gene_interval1 and gene_interval2:
				if gene_interval1.chr == gene_interval2.chr:
					self.gene_interval_distance = max(gene_interval1.start, gene_interval2.start) - min(gene_interval1.end, gene_interval2.end)

			#switch pos order if necessary
			if switch_pos:
				self.chr1, self.chr2 = self.chr2, self.chr1
				self.pos1, self.pos2 = self.pos2, self.pos1
				self.strand1, self.strand2 = self.strand2, self.strand1
				self.t_gene1, self.t_gene2 = self.t_gene2, self.t_gene1
				self.t_area1, self.t_area2 = self.t_area2, self.t_area1


		#print max_t1[0], max_t2[0]

		'''
		gene_order.append(",".join(list1) + ">>" + ",".join(list2))
		gene_order.append(",".join(type1) + ">>" + ",".join(type2))
		gene_order.append(",".join(id1) + ">>" + ",".join(id2))
		gene_order.append(category)
		'''
		#gene1 = "NA" if not self.bpann1 else self.bpann1.gene_name
		#gene2 = "NA" if not self.bpann2 else self.bpann2.gene_name
		return ""
	# based on given gene annotations re-annotate cff fusions, infer possible up/downstream genens, try to fill in strand if info missing
	def ann_gene_order(self, gene_ann):
		gene_order = []
		
		matched_genes1 = gene_ann.map_pos_to_genes(self.chr1, self.pos1)
		matched_genes2 = gene_ann.map_pos_to_genes(self.chr2, self.pos2)
		
		a = {} # forward strand gene at pos1
		c = {} # backward strand gene at pos1
		b = {} # forward strand gene at pos2
		d = {} # backward strand gene at pos2
		for gene in matched_genes1:
			if gene.strand == "f":
				a.setdefault(gene.gene_name, []).append(gene)
			else:
				c.setdefault(gene.gene_name, []).append(gene)
		for gene in matched_genes2:
			if gene.strand == "f":
				b.setdefault(gene.gene_name, []).append(gene)
			else:
				d.setdefault(gene.gene_name, []).append(gene)
		# for tools do not provide defuse-style strand, regenerate strands, this is assuming that gene1 = 5 prime gene and gene2 = 3 prime gene
		if self.strand1 == "NA" or self.strand2 == "NA":
			gene_interval1 = ""
			gene_interval2 = ""
			for sep in [",", "/"]:
				for gene_name in (self.t_gene1).split(sep):
					if not gene_interval1:
						gene_interval1 = gene_ann.get_gene_interval(gene_name)
					else:
						break
				for gene_name in (self.t_gene2).split(sep):
					if not gene_interval2:
						gene_interval2 = gene_ann.get_gene_interval(gene_name)
					else:
						break
				
			# Gene Strand : Defuse Strand
			# + >> + : + -
			# + >> - : + +
			# - >> + : - -
			# - >> - : - +
			if gene_interval1 and gene_interval2:
				if gene_interval1.strand == "f":
					self.strand1 = "+"
				else:
					self.strand1 = "-"
				if gene_interval2.strand == "f":
					self.strand2 = "-"
				else:
					self.strand2 = "+"
			else:
				# failed to fill in strand, return list with warnning info
				#gene_order.append("Strand_filling_failed")
				return gene_order

				
		# gene_order includes: 5' gene >> 3' gene, 5' gene type >> 3' gene type, 5' coding gene idx >> 3' coding gene inx, category
		if self.strand1 == "+" and self.strand2 == "+":
			gene_order = self.__check_gene_pairs(a, d, gene_ann, False)
			gene_order += self.__check_gene_pairs(b, c, gene_ann, True)
		elif self.strand1 == "+" and self.strand2 == "-":
			gene_order = self.__check_gene_pairs(a, b, gene_ann, False)
			gene_order += self.__check_gene_pairs(d, c, gene_ann, True)
		elif self.strand1 == "-" and self.strand2 == "+":
			gene_order = self.__check_gene_pairs(c, d, gene_ann, False)
			gene_order += self.__check_gene_pairs(b, a, gene_ann, True)
		elif self.strand1 == "-" and self.strand2 == "-":
			gene_order = self.__check_gene_pairs(c, b, gene_ann, False)
			gene_order += self.__check_gene_pairs(d, a, gene_ann, True)
		#self.reann_gene_order1, self.reann_gene_type1, self.reann_gene_index1, self.reann_category1, self.reann_gene_order2, self.reann_gene_type2, self.reann_gene_index2, self.reann_category2 = gene_order
		
		'''
		# check whether pos on mapped genes boundaries
		on_boundary1 = False
		close_to_boundary1 = False
		on_boundary2 = False
		close_to_boundary2 = False
		for gene in matched_genes1:
			#if not gene.type == "cds":
			#	continue
			if gene.is_on_boundary:
				on_boundary1 = True
				close_to_boundary1 = True
				break
			if gene.close_to_boundary:
				close_to_boundary1 = True
		for gene in matched_genes2:
			#if not gene.type == "cds":
			#	continue
			if gene.is_on_boundary:
				on_boundary2 = True
				close_to_boundary2 = True
				break
			if gene.close_to_boundary:
				close_to_boundary2 = True
		self.boundary_info = "\t".join(map(str, [on_boundary1, close_to_boundary1, on_boundary2, close_to_boundary2]))
		gene_order += [on_boundary1, close_to_boundary1, on_boundary2, close_to_boundary2]
		return gene_order
		'''

		'''
		# test pre/next ann
		print "pre:"
		print "cur", self.bpann1.tostring()	
		for ann in gene_ann.get_previous_ann(self.bpann1):
			print ann.tostring()
		print "next:"
		print "cur", self.bpann2.tostring()	
		for ann in gene_ann.get_next_ann(self.bpann2):
			print ann.tostring()
		'''	

	# realign breakpoints of this fusion to the left most, not finished, how to define "left" when genes are on different chrs 
	def left_aln_fusion_bp(self, refs):
		# provided reference file lacks fusion chr
		if not (self.chr1 in refs.references and self.chr2 in refs.references):
			return (-1, -1)
		rlen = 10
		#refs = pysam.FastaFile(ref_file)
		
		#pysam use 0-based coordinates, cff use 1-based coordinates
		if self.strand1 == "+" and self.strand2 == "+":
			#up_seq = sequtils.rc_seq(refs.fetch(self.chr1, self.pos1-rlen, self.pos1+rlen), "r")
			#down_seq = sequtils.rc_seq(refs.fetch(self.chr2, self.pos2-rlen, self.pos2+rlen), "c")
			up_seq = refs.fetch(self.chr1, self.pos1-rlen, self.pos1+rlen)
			down_seq = sequtils.rc_seq(refs.fetch(self.chr2, self.pos2-rlen, self.pos2+rlen), "rc")
		elif self.strand1 == "+" and  self.strand2 == "-":
			#up_seq = sequtils.rc_seq(refs.fetch(self.chr1, self.pos1-rlen, self.pos1+rlen), "r")
			#down_seq = sequtils.rc_seq(refs.fetch(self.chr2, self.pos2-rlen-1, self.pos2-1+rlen), "r")
			up_seq = refs.fetch(self.chr1, self.pos1-rlen, self.pos1+rlen)
			down_seq = refs.fetch(self.chr2, self.pos2-rlen-1, self.pos2-1+rlen)
		elif self.strand1 == "-" and  self.strand2 == "+":
			#up_seq = refs.fetch(self.chr1, self.pos1-1-rlen, self.pos1+rlen-1)
			#down_seq = refs.fetch(self.chr2, self.pos2-rlen, self.pos2+rlen)
			down_seq = refs.fetch(self.chr1, self.pos1-1-rlen, self.pos1+rlen-1)
			up_seq = refs.fetch(self.chr2, self.pos2-rlen, self.pos2+rlen)
		elif self.strand1 == "-" and  self.strand2 == "-":
			#up_seq = sequtils.rc_seq(refs.fetch(self.chr1, self.pos1-1-rlen, self.pos1+rlen-1), "c")
			#down_seq = sequtils.rc_seq(refs.fetch(self.chr2, self.pos2-rlen-1, self.pos2-1+rlen), "r")
			#up_seq = sequtils.rc_seq(refs.fetch(self.chr1, self.pos1-rlen, self.pos1+rlen), "rc")
			#down_seq = refs.fetch(self.chr2, self.pos2-rlen-1, self.pos2-1+rlen)
			down_seq = refs.fetch(self.chr1, self.pos1-rlen-1, self.pos1+rlen-1)
			up_seq = sequtils.rc_seq(refs.fetch(self.chr2, self.pos2-rlen-1, self.pos2+rlen-1), "rc")
		else:
			print >> sys.stderr, "Unknown strand:", self.strand1, self.strand2
			sys.exit(1)
		
		if len(up_seq) < rlen or len(down_seq) < rlen:
			print >> sys.stderr, "Warnning: reference sequence cannot be fetched."
			print >> sys.stderr, self.tostring()
		
			return (-1, -1)

		i = rlen - 1
		while i >= 0:
			if up_seq[i].upper() == down_seq[i].upper():
				i -= 1
			else:
				break
		j = rlen
		while j < 2*rlen:
			if up_seq[j].upper() == down_seq[j].upper():
				j += 1
			else:
				break

		print up_seq.upper()[0:rlen], up_seq.lower()[rlen:]
		print down_seq.lower()[0:rlen], down_seq.upper()[rlen:]
		return (rlen-1-i, j-rlen)
	# aln fusion breakpoint to exon/utr boundary
	def aln_seq_to_boundary(self, seq1, splice_seq1, splice_seq2, seq2, distance_to_bnd, strand):
		'''
		print "seq1:", seq1[-10:]
		print "splice seq1:", splice_seq1
		print "splice seq2:", splice_seq2
		print "seq2:", seq2[:10]
		'''
		step = 0 
		gt_ag_step = -1
		if distance_to_bnd > 0: # move right
			for c in seq2:
				if step == len(splice_seq1):
					break
				splice_site1 = splice_seq1[step:step+2]
				splice_site2 = splice_seq2[len(splice_seq2)-2+step:] + seq2[max(0, step-2):step]
				if splice_site1 == "GT" and splice_site2 == "AG":
					gt_ag_step = step
				if step == distance_to_bnd:
					break
				if c == splice_seq1[step]:
					step += 1
				else:
					break

		else: # move left
			for c in seq1[::-1]:
				if step == len(splice_seq2):
					break
				splice_site1 = seq1[len(seq1)-step:len(seq1)-step+2] + splice_seq1[:max(0, 2-step)]
				splice_site2 = splice_seq2[-2-step:len(splice_seq2)-step]
				if splice_site1 == "GT" and splice_site2 == "AG":
					gt_ag_step = step
				if step == abs(distance_to_bnd):
					break
				if c == splice_seq2[-1-step]:
					step += 1
				else:
					break
		'''
		print distance_to_bnd
		print step
		print gt_ag_step
		'''
		if step == abs(distance_to_bnd): # shift to boundary
			self.splice_site1 = splice_site1
			self.splice_site2 = splice_site2
			if strand == "f":
				self.pos1 += distance_to_bnd
				self.pos2 += distance_to_bnd
			else:
			
				self.pos1 -= distance_to_bnd
				self.pos2 -= distance_to_bnd
		elif gt_ag_step > 0: #shift to GT-AG breakpoint
			step = gt_ag_step
			self.splice_site1 = splice_seq1[step:step+2] 
			self.splice_site2 = splice_seq2[len(splice_seq2)-2+step:] + seq2[max(0, step-2):step]
			if distance_to_bnd < 0: gt_ag_step = -gt_ag_step
			if strand == "f":
				self.pos1 += gt_ag_step
				self.pos2 += gt_ag_step
			else:
			
				self.pos1 -= gt_ag_step
				self.pos2 -= gt_ag_step
		else: # do not shift
			step = 0
			self.splice_site1 = splice_seq1[step:step+2] 
			self.splice_site2 = splice_seq2[len(splice_seq2)-2+step:] + seq2[max(0, step-2):step]
			
			
		#print self.splice_site1, self.splice_site2
	# get coding sequence up/downstream of fusion breakpoint. trans_ann_list is a list of GeneBed for the same transcript. ref is an opened pysam fasta
	def get_transcript_seq(self, trans_ann_list, bp, ref, gene_order):
		if not trans_ann_list:
			print "empty trans_ann_list"
			return "", "", -1, -100, "", ""
		#else:
		#	print "ann num:", len(trans_ann_list)
		'''
		if coding:
			ann_list = sorted(filter(lambda bpann:bpann.type=="cds", trans_ann_list), key = lambda x:x.start)
		else:
			ann_list = sorted(filter(lambda bpann:bpann.type!="intron", trans_ann_list), key = lambda x:x.start)
		'''	
		seqs = []
		is_fw = (trans_ann_list[0].strand == "f")
		is_head = (gene_order == "head")
		splice_site = "NA"
		trans_id = trans_ann_list[0].transcript_id

		bp_in_coding = False # breakpoint of head gene in coding exon

		#get last exon idx for current transcript
		exon_ann_list = filter(lambda bpann:bpann.type=="cds", trans_ann_list)
		if exon_ann_list:
			if is_fw:
				last_exon_idx = exon_ann_list[-1].idx
			else:
				last_exon_idx = exon_ann_list[0].idx
		else:
			last_exon_idx = -1
		
		last_exon_end = -1 # for tail seq, get the position of last basepair of last exon in seq2
		distance_to_bnd = -100 #distance to exon/utr boundary, if breakpoint in exon/utr it equals to bnd - bp, thus positive, otherwise negative
		
		if is_head:
			if is_fw: # head gene on forward strand
				for ann in trans_ann_list:
					# breakpoint in current bed annotation
					if ann.start <= bp <= ann.end:
						seq = ref.fetch(format_chr(ann.chr, ref), ann.start-1, bp)
						seqs.append(seq) 
						splice_site = ref.fetch(format_chr(ann.chr, ref), bp, bp+20)
						if ann.type == "cds" or len(seq) < 10: # breakpoint in cds, or close to cds
							bp_in_coding = True
						if ann.type == "intron":
							distance_to_bnd = ann.start - bp - 1
						else:
							distance_to_bnd = ann.end - bp
						bp_ann = ann
						
						break
					elif ann.type == "cds":
						seqs.append(ref.fetch(format_chr(ann.chr, ref), ann.start-1, ann.end))
			else:
				for ann in trans_ann_list:
					# breakpoint in current bed annotation
					if ann.start <= bp <= ann.end:
						seq = ref.fetch(format_chr(ann.chr, ref), bp-1, ann.end)
						seqs.append(seq) 
						
						splice_site = ref.fetch(format_chr(ann.chr, ref), bp-21, bp-1)
						if ann.type == "cds" or len(seq) < 10: # breakpoint in cds, or close to cds
							bp_in_coding = True
						if ann.type == "intron":
							distance_to_bnd = bp - ann.end - 1
						else:
							distance_to_bnd = bp - ann.start
						bp_ann = ann
					elif ann.start > bp and ann.type == "cds": # for now do not consider intonic breakpoint
						seqs.append(ref.fetch(format_chr(ann.chr, ref), ann.start-1, ann.end))

		else:
			if is_fw: # tail gene on forward strand
				for ann in trans_ann_list:
					if ann.type=="cds" and ann.idx == last_exon_idx:
						last_exon_end = len("".join(seqs)) + (ann.end  - ann.start + 1)
					# breakpoint in current bed annotation
					if ann.start <= bp <= ann.end:
						splice_site = ref.fetch(format_chr(ann.chr, ref), bp-21, bp-1)
						seqs.append(ref.fetch(format_chr(ann.chr, ref), bp-1, ann.end)) 
						if ann.type == "cds":
							bp_in_coding = True
						bp_ann = ann
					elif ann.start > bp and ann.type != "intron": # for now do not consider intonic breakpoint
						seqs.append(ref.fetch(format_chr(ann.chr, ref), ann.start-1, ann.end))
			else:
				for ann in trans_ann_list:
					# breakpoint in current bed annotation
					if ann.start <= bp <= ann.end:
						seqs.append(ref.fetch(format_chr(ann.chr, ref), ann.start-1, bp)) 
						splice_site = ref.fetch(format_chr(ann.chr, ref), bp, bp+20)
						if ann.type == "cds":
							bp_in_coding = True
						bp_ann = ann
						break
					elif ann.type != "intron": # for now do not consider intonic breakpoint
						seqs.append(ref.fetch(format_chr(ann.chr, ref), ann.start-1, ann.end))

		'''
		# head gene on forward strand, all exons before bp and for the ann contains bp, output seq from start to bp
		if (is_head and is_fw) or (not is_head and not is_fw) :
			for ann in trans_ann_list:
				if ann.start <= bp <= ann.end:
					seqs.append(ref.fetch(ann.chr, ann.start-1, bp)) 
					splice_site = ref.fetch(ann.chr, bp, bp+2)
					if ann.type == "cds":
						bp_in_coding = True
					if ann.idx == last_exon_idx:
						
					break
				elif ann.type == "cds":
					seqs.append(ref.fetch(ann.chr, ann.start-1, ann.end))
		else: # head gene on reverse strand, start from bp to the end of ann that contains it, then all exons after bp
			for ann in trans_ann_list:
				if ann.start <= bp <= ann.end:
					splice_site = ref.fetch(ann.chr, bp-3, bp-1)
					seqs.append(ref.fetch(ann.chr, bp-1, ann.end)) 
					if ann.type == "cds":
						bp_in_coding = True
				elif ann.start > bp and ann.type != "intron": # for now do not consider intonic breakpoint
					seqs.append(ref.fetch(ann.chr, ann.start-1, ann.end))
		if is_head and not bp_in_coding:		
			print "head gene breakpoint not in cds", bp, trans_id
			return "", "", -1, -100
		'''
		#trans_seq = "___".join(seqs)
		trans_seq = "".join(seqs).upper()
		splice_site = splice_site.upper()		
		if not is_fw:
			trans_seq = sequtils.rc_seq(trans_seq, "rc")
			splice_site = sequtils.rc_seq(splice_site, "rc")
		'''
		if not is_head and last_exon_end > 0:
			print "last exon end:", last_exon_end
			print "last codon:", trans_seq[last_exon_end-3:last_exon_end]
		'''
		return trans_seq, splice_site, last_exon_end, distance_to_bnd, bp_ann, bp_in_coding	
			


	# for given cff fusion, test all combination of head/tail gene sequences to check their codons. Can be used to classiy readthrough more strictly by requiring the fusion seq in-frame and reach the end of tail gene		
	def check_codon(self, gene_ann, ref_file):
		head_seqs = []
		tail_seqs = []
		ref = pysam.FastaFile(ref_file)	
		# get all transcripts for head gene	
		head_gene_intervals = gene_ann.get_gene_interval(self.reann_gene1)
		tids = head_gene_intervals.transcript_ids
		if not head_gene_intervals.is_coding:
			print >> sys.stderr, "head gene is not coding"
			return
			
		if not tids:
			print >> sys.stderr, "Head gene not in annotation:", self.reann_gene1
			return

		for transcript_id in tids:
			bpann_list = gene_ann.get_transcripts_ann(transcript_id) # all GeneBed annotations of transcript_id
			if not bpann_list[0].start <= self.pos1 <= bpann_list[-1].end: #skip the transcript that doesnt include the breakpoint
				continue
			seq, splice_site, last_exon_end, distance_to_bnd, bp_ann, bp_in_coding =  self.get_transcript_seq(bpann_list, self.pos1, ref, "head")
			head_seqs.append((seq, splice_site, last_exon_end, transcript_id, distance_to_bnd, bp_ann, bp_in_coding))
		# get all transcripts for tail gene	
		tids = (gene_ann.get_gene_interval(self.reann_gene2)).transcript_ids
		if not tids:
			# no tail gene, make a fake bed annotation, extract 1k sequence downstream the breakpint
			fake_ann = GeneBed("")
			fake_ann_len = 1000
			fake_ann.chr = self.chr2
			fake_ann.type = "intergenic"
			fake_ann.idx = 0
			fake_ann.start = self.pos2 - fake_ann_len
			fake_ann.end = self.pos2 + fake_ann_len

			# according to head gene strand, infer fake_ann strand
			if self.strand1 == self.strand2: # head gene and tail gene on different strand	 
				if head_gene_intervals.strand == "f":
					fake_ann.strand = "r"
				else:
					fake_ann.strand = "f"
			else: # head and tail gene on same strand
				fake_ann.strand = head_gene_intervals.strand
			seq, splice_site, last_exon_end, distance_to_bnd, bp_ann, bp_in_coding = self.get_transcript_seq([fake_ann], self.pos2, ref, "tail")
			tail_seqs.append((seq, splice_site, last_exon_end, fake_ann.transcript_id, distance_to_bnd, bp_ann, bp_in_coding))
		for transcript_id in tids:
			bpann_list = gene_ann.get_transcripts_ann(transcript_id) # all GeneBed annotations fo transcript_id
			if not bpann_list[0].start <= self.pos2 <= bpann_list[-1].end:
				continue
			seq, splice_site, last_exon_end, distance_to_bnd, bp_ann, bp_in_coding = self.get_transcript_seq(bpann_list, self.pos2, ref, "tail")
			tail_seqs.append((seq, splice_site, last_exon_end, transcript_id, distance_to_bnd, bp_ann, bp_in_coding))
		# search first stop codon in tail gene sequence, stop codon: TAA, TAG, TGA
		stop_codons= set(["TAA", "TAG", "TGA"])
		seq_set = set()
		# among all the combinations of head and tail sequences, 1. check if the fusion can reach any transcript's end (first stop codon is the tail transcript's stop codon); 2. if not 1, return the smallest stop condon position.
		infered_fusion_seq_info = ("", False, -1, "", "", -100, "", 0) # infered seq, reach last exon, first stop codon pos, upstream transcript id, downstream gene transcript_id, distance_to_bnd, score
		for seq1, splice_site1, last_exon_end1, transcript_id1, distance_to_bnd1, bp_ann1, bp_in_coding1 in head_seqs:
			if not seq1:
				continue
			shift1 = len(seq1) % 3
			for seq2, splice_site2, last_exon_end2, transcript_id2, distance_to_bnd2, bp_ann2, bp_in_coding2 in tail_seqs:
				i = 0
				shift_seq2 = seq2[(3-shift1)%3:]
				# convert shift_seq2 to triplets
				shift_seq2_triplets = re.findall(r'.{3}', shift_seq2)
				# get uniq triplets set
				shift_seq2_codon_set = set(shift_seq2_triplets)
				# compare to stop codons
				shift_seq2_stop_conons =  stop_codons & shift_seq2_codon_set
				# stop codon found, get the first one
				if shift_seq2_stop_conons:
					first_stop_idx = len(shift_seq2_triplets)
					first_stop_codon = "NA"
					for codon in shift_seq2_stop_conons:
						stop_idx = shift_seq2_triplets.index(codon)	
						first_stop_idx = stop_idx
						first_stop_codon = codon
						#if stop_idx < first_stop_idx:
						#	first_stop_idx = stop_idx
						#	first_stop_codon = codon
					tail_gene_seq = seq2[:(3-shift1)%3 + (first_stop_idx + 1)*3]	
					#if not tail_gene_seq in seq_set:
					update = False
					score = 0 # current 
					is_inframe = False
					if bp_in_coding1 and bp_in_coding2 and (len(seq1) + len(seq2)) % 3 == 0: # head transcript is coding and inframe
						is_inframe = True
						score += 10000

					if abs(infered_fusion_seq_info[5]) > abs(distance_to_bnd1): # current combination closer to annotation boundary
						score += 1000
					if bp_ann1.type == "cds": # head transcript invole exon
						score += 500
					if bp_ann2.type == "cds": 
						score += 300
					if bp_ann1.type == "utr5" or  bp_ann1.type == "utr3": # head transcript involve utr
						score += 200
					if bp_ann2.type == "utr5" or bp_ann2.type == "utr3": 
						score += 100
					'''
					if abs(infered_fusion_seq_info[5]) < abs(distance_to_bnd1):
						score1 += 1000
					if abs(infered_fusion_seq_info[5]) > abs(distance_to_bnd1):
						score2 += 1000
					
					if bp_ann1.type == "cds": # saved combination reaches last exon
						score1 += 500
					if bp_ann2.type == "cds": # current combination reaches last exon:
						score2 += 400

					if bp_ann1.type == "utr5" or  bp_ann1.type == "utr3": # saved combination reaches last exon
						score1 += 200
					if bp_ann2.type == "utr5" or bp_ann2.type == "utr3": # current combination reaches last exon:
						score2 += 100
					'''
					if gene_ann.transcript_is_coding(transcript_id2): # current transcript has cds
						score += 10
					if first_stop_idx <= int(infered_fusion_seq_info[2]): # saved combination has a smaller stop codon position
						score += 1
					update = True if score > int(infered_fusion_seq_info[7]) else False
					if update:
						#print "Scores:", score, infered_fusion_seq_info[5], distance_to_bnd1, len(seq1), len(seq2), infered_fusion_seq_info[6]		
						#print "stop codon found:", first_stop_codon, first_stop_idx
						#print "length:", len(seq1), len(seq2), len(seq2[:(3-shift1)%3 + (first_stop_idx + 1)*3]), last_exon_end2
				
						# infered fusion sequence, ends at first stop codon
						infered_fusion_seq = seq1 + "_" + splice_site1 + "__" + splice_site2 + "_" + seq2[:(3-shift1)%3 + (first_stop_idx + 1)*3]
						#infered_fusion_seq_info = (infered_fusion_seq, len(seq2[:(3-shift1)%3 + (first_stop_idx + 1)*3])==last_exon_end2, first_stop_idx, transcript_id1, transcript_id2, distance_to_bnd1, bp_ann1, score)
						infered_fusion_seq_info = (infered_fusion_seq, is_inframe, first_stop_idx, transcript_id1, transcript_id2, distance_to_bnd1, bp_ann1, score)
					#seq_set.add(tail_gene_seq)	
					#print "trans id:", transcript_id2
					#print "stop codon found:", first_stop_codon, first_stop_idx
					#print "length:", len(seq1), len(seq2[:(3-shift1)%3 + (first_stop_idx + 1)*3]), last_exon_end2
					#print "seq1:", seq1 + "|" + splice_site1
					#print "seq2:", splice_site2 + "|" + seq2[:(3-shift1)%3 + (first_stop_idx + 1)*3]
					
					#seq_set.add((seq1, splice_site1, seq2, splice_site2, first_stop_idx, first_stop_codon))
					
				else:
					#print "no stop codon found:", transcript_id2
					#print [ann.type for ann in gene_ann.get_transcripts_ann(transcript_id2)]
					pass
					
		#print infered_fusion_seq_info
		is_inframe = infered_fusion_seq_info[1]
		seqs = infered_fusion_seq_info[0].split("_")
		distance_to_bnd = infered_fusion_seq_info[5]
		bp_ann = infered_fusion_seq_info[6]
		if len(seqs) == 5:
			#print "before:", self.pos1, self.pos2
			#print infered_fusion_seq_info[0]
			self.aln_seq_to_boundary(seqs[0], seqs[1], seqs[3], seqs[4], distance_to_bnd, bp_ann.strand)

			#print "after:", self.pos1, self.pos2
			self.is_inframe = is_inframe
			#print seqs[0][-10:]
			#print seqs[4]
			#print len(seqs[4])
		#for seq in seqs: print len(seq),
		#print
		
		#print len(head_seq), len(tail_seq), len(head_seq) + len(tail_seq)
		'''
		for seq1, splice_site1, seq2, splice_site2, first_stop_idx, first_stop_codon in seq_set:
			continue
			shift1 = len(seq1) % 3
			print "stop codon found:", first_stop_codon, first_stop_idx
			print "length:", len(seq1), len(seq2[:(3-shift1)%3 + (first_stop_idx + 1)*3])
			print "seq1:", seq1 + "|" + splice_site1
			print "seq2:", splice_site2 + "|" + seq2[:(3-shift1)%3 + (first_stop_idx + 1)*3]
			#print "seq1:", ",".join(re.findall(r'.{3}', seq1)) + "," + seq1[-(len(seq1) % 3):] + "|" + splice_site1
			#print "seq2:", splice_site2 + "|" + seq2[:(3-shift1)%3] + "," + ",".join(shift_seq2_triplets[:first_stop_idx+1])
		'''	
		ref.close()
		
		
class GeneIntervals():
	def __init__(self, bed_ann_list):
		'''
		self.gene_name = gene_name
		self.chr = chr
		self.start = start
		self.end = end
		self.strand = strand
		self.is_coding = is_coding
		'''
		self.load(bed_ann_list)
	def check_bed_ann_list(self, bed_ann_list):
		flag = True
		if not bed_ann_list:
			flag = False
		else:
			if len(set([(a.gene_name, a.chr, a.strand) for a in bed_ann_list])) > 1:
				flag = False
			else:
				max_start = max([a.start for a in bed_ann_list])
				min_end = min([a.end for a in bed_ann_list])
				# Gene has two annotations more than 1Mb away from each other
				if max_start - min_end > 5000000:
					flag = False
			if not flag:
				print >> sys.stderr, "Warnning: Input gene annotations include multiple chr, strand, or regions (5Mb away). Skipping current gene annotation."
				print >> sys.stderr, set([(a.gene_name, a.chr, a.strand) for a in bed_ann_list])

		return flag
	def load(self, bed_ann_list):
		if not self.check_bed_ann_list(bed_ann_list):
			self.gene_name = "NA"
			self.chr = "NA"
			self.start = -1
			self.end = -1
			self.strand = "NA"
			self.is_coding = False
			self.is_contradictory = True # contradictory gene annotation, will not be used
			self.transcript_ids = []
		else:
			self.gene_name = bed_ann_list[0].gene_name
			self.chr = bed_ann_list[0].chr
			self.strand= bed_ann_list[0].strand
			self.start = min([a.start for a in bed_ann_list])
			self.end = max([a.end for a in bed_ann_list])
			self.is_coding = True if "cds" in [a.type for a in bed_ann_list] else False
			self.is_contradictory = False
			self.transcript_ids = list(set([genebed.transcript_id for genebed in bed_ann_list]))
		'''
		if not bed_ann_list:
			self.gene_name = "NA"
			self.chr = "NA"
			self.start = -1
			self.end = -1
			self.strand = "NA"
			self.is_coding = False
			#print >> sys.stderr, "Empty GeneBed annotation."
			#sys.exit(1)
		else:
			# warnning when same gene on different chr/strand, load  the first annotation
			if len(set([(a.gene_name, a.chr, a.strand) for a in bed_ann_list])) > 1:
				print >> sys.stderr, "Warnning: GeneBed annotation includes multiple genes. Only the first annotation will be used."
				print >> sys.stderr, set([(a.gene_name, a.chr, a.strand) for a in bed_ann_list])
				#sys.exit(1)
			self.gene_name = bed_ann_list[0].gene_name
			self.chr = bed_ann_list[0].chr
			self.strand= bed_ann_list[0].strand

			self.start = min([a.start for a in filter(lambda x:x.chr==self.chr, bed_ann_list)])
			self.end = max([a.end for a in filter(lambda x:x.chr==self.chr, bed_ann_list)])
			#self.start = min([a.start for a in bed_ann_list])
			#self.end = max([a.end for a in bed_ann_list])
			#self.is_coding = True if "cds" in [a.type for a in bed_ann_list] else False
			self.is_coding = True if "cds" in [a.type for a in filter(lambda x:x.chr==self.chr, bed_ann_list)] else False
		'''		
	
	def overlap(self, interval2):
		if self.chr == interval2.chr and min(self.end, interval2.end) - max(self.start, interval2.start) > 0:
			return True
		else:
			return False
	def merge(self, interval2):
		if self.chr != interval2.chr or self.strand != interval2.strand:
			print >> sys.stderr, "Warning: intervals are on different chr/strand."
			print >> sys.stderr, self.gene_name, interval2.gene_name
			return self
		else:
			#new_interval = GeneIntervals("Merged_" + self.gene_name + "_" + interval2.gene_name, self.chr, min(self,start, interval2.start), max(self.end, interval2.end), self.strand, False)
			new_interval = copy.copy(interval2)
			new_interval.gene_name = "Merged_" + self.gene_name + "_" + interval2.gene_name
			new_interval.chr =  self.chr
			new_interval.start =  min(self.start, interval2.start)
			new_interval.end =  max(self.end, interval2.end)
			new_interval.strand = self.strand
			new_interval.transcript_ids = list(set(self.transcript_ids + interval2.transcript_ids))
			# for merged gene intervals is_coding has no sense
			new_interval.is_coding = False
			
			return new_interval			
# Deprecated, use GeneBed instead
class BreakpointAnnotation:
	#tuple = (chr, start, end, transcript, type, idx, strand, gene_name)
	def __init__(self, chr, start, end, transcript_name, type, idx, strand, gene_name):
		self.__attrs = []
		self.chr = chr 			
		self.start = start
		self.end = end
		self.transcript_name = transcript_name
		self.type = type # utr/cds/intron
		self.idx = idx
		self.strand = strand
		self.gene_name = gene_name
		self.is_on_boundary = False
		self.__attrs.append(chr)
		self.__attrs.append(str(start))
		self.__attrs.append(str(end))
		self.__attrs.append(transcript_name)
		self.__attrs.append(type)
		self.__attrs.append(str(idx))
		self.__attrs.append(strand)
		self.__attrs.append(gene_name)
	def tostring(self):
		return  "\t".join(self.__attrs)
#Load bed format gene annotation, current support knowngene.bed's format, map given genomic loactions to genens, return matched gene list
class GeneAnnotation():
	__gene_starts = {}
	__genes = {}
	__max_diff = 1000000
	__gene_intervals = {}
	__transcripts_ann = {}
	__gene_name_id_map = {}
	__gene_name_idx_map = {}
	__coding_gene_list = []
	# gene_ann_bed
	# chr1    24683494	24685032	ENST00000003583 utr3    0       r       STPG1   ENSG00000001460
	def __init__(self, gene_ann_bed):
		if gene_ann_bed != "":
			self.load_gene_bed(gene_ann_bed)
			self.load_gene_intervals(gene_ann_bed)
			self.build_coding_gene_list()
			'''
			print >> sys.stderr, "mem of gene_starts:", sys.getsizeof(self.__gene_starts)
			print >> sys.stderr, "length:", len(self.__gene_starts.keys())
			print >> sys.stderr, "unit size:", sys.getsizeof(self.__gene_starts[self.__gene_starts.keys()[0]])
				
			print >> sys.stderr, "mem of genes:", sys.getsizeof(self.__genes)
			print >> sys.stderr, "length:", len(self.__genes.keys())
			print >> sys.stderr, "unit size:", sys.getsizeof(self.__genes[self.__genes.keys()[0]])

			print >> sys.stderr, "mem of gene_intervals:", sys.getsizeof(self.__gene_intervals)
			print >> sys.stderr, "length:", len(self.__gene_intervals.keys())
			print >> sys.stderr, "unit size:", sys.getsizeof(self.__gene_intervals[self.__gene_intervals.keys()[0]].start)
			'''
			self.load_transcripts_ann(gene_ann_bed)
	# for each transcript save all its cds,intron, utr annotations in dictionary __transcripts_ann
	def load_transcripts_ann(self, gene_ann_bed):
		start_time = time.time()
		for line in open(gene_ann_bed, "r"):
			ann = GeneBed(line)
			key = ann.transcript_id
			self.__transcripts_ann.setdefault(key, []).append(ann)
		print >> sys.stderr, "Transcript annotations loaded."
		print >> sys.stderr, time.time() - start_time, "sec. elapsed."	
	
	def get_transcripts_ann(self, trans_id):
		if trans_id in self.__transcripts_ann:
			return self.__transcripts_ann[trans_id]
		else:
			return []	
	def transcript_is_coding(self, trans_id):
		if trans_id in self.__transcripts_ann:
			if filter(lambda ann:ann.type=="cds" ,self.__transcripts_ann[trans_id]):
				return True
		return False	

	def get_previous_ann(self, bpann):
		bpann_list = self.get_transcripts_ann(bpann.transcript_id)
		pre_bpann_list = []
		for ann in bpann_list:
			if bpann.strand == "f":
				if ann.idx == bpann.idx - 1:
					pre_bpann_list.append(ann)
			elif bpann.strand == "r":
				if ann.idx == bpann.idx + 1:
					pre_bpann_list.append(ann)
		return pre_bpann_list	
	def get_next_ann(self, bpann):
		bpann_list = self.get_transcripts_ann(bpann.transcript_id)
		next_bpann_list = []
		for ann in bpann_list:
			if bpann.strand == "f":
				if ann.idx == bpann.idx + 1:
					next_bpann_list.append(ann)
			elif bpann.strand == "r":
				if ann.idx == bpann.idx - 1:
					next_bpann_list.append(ann)
		return next_bpann_list	
	def get_gene_id(self, gene):
		if gene in self.__gene_name_id_map:
			return self.__gene_name_id_map[gene]
		else:
			return "NA"
	# for each gene use minimal start and max end of all its transcripts as its interval, same intervals of all genes in an dictionary __gnen_intervals
	def load_gene_intervals(self, gene_ann_bed):
		start_time = time.time()
		n1 = 0
		n2 = 0
		tmp_dict = {}
		# use gene name as key to bulid a dict for gene_ann_bed, all annotations for the same gene are saved in a list which will be used for getting intervals
		for line in open(gene_ann_bed, "r"):
			bed_ann = GeneBed(line)
			if self.filter_gene_ann(bed_ann):
				continue
			tmp_dict.setdefault(bed_ann.gene_name, []).append(bed_ann)
		for gene_name in tmp_dict:
			self.__gene_intervals.setdefault(gene_name, GeneIntervals(tmp_dict[gene_name]))
		print >> sys.stderr, "Gene intervals loaded."
		print >> sys.stderr, time.time() - start_time, "sec. elapsed."	
		#print >> sys.stderr, "#gene on diff chr:", n1, "#gene on diff strand:", n2

	# index foward/backward genes to decide whether two genes are next to each other (possible readthrough)
	def build_coding_gene_list(self):
		for gene_name in self.__gene_intervals:
			interval = self.__gene_intervals[gene_name]
			if interval.is_coding:
				self.__coding_gene_list.append(interval)
		self.__coding_gene_list.sort(key = lambda i:(i.chr, i.start))

		i_f = 0 # idx of forward strand gene
		i_r = 0 # idx of reverse strand gene
		#pre_interval_f = GeneIntervals("None", "chr0", 0, 0, "+", False)
		#pre_interval_r = GeneIntervals("None", "chr0", 0, 0, "-", False)
		pre_interval_f = GeneIntervals([])
		pre_interval_r = GeneIntervals([])
		for interval in self.__coding_gene_list:
			if interval.strand == "+" or interval.strand == "f":
				if interval.overlap(pre_interval_f):
					pre_interval_f = interval.merge(pre_interval_f)
				else:
					pre_interval_f = interval
					i_f += 1
				self.__gene_name_idx_map.setdefault(interval.gene_name, str(i_f) + "_" + interval.strand)
			elif interval.strand == "-" or interval.strand == "r":
				if interval.overlap(pre_interval_r):
					pre_interval_r = interval.merge(pre_interval_r)
				else:
					pre_interval_r = interval
					i_r += 1
				self.__gene_name_idx_map.setdefault(interval.gene_name, str(i_r) + "_" + interval.strand)
	def print_coding_gene_idx(self):
		for key in self.__gene_name_idx_map:
			print key, self.__gene_name_idx_map[key]
	# whether a gene include coding exon (cds)
	def is_coding(self, gene_name):
		return self.__gene_intervals[gene_name].is_coding 
	# where a gene has contradictory annotations
	def is_contradictory(self, gene_name):
		return self.__gene_intervals[gene_name].is_contradictory
	def get_coding_gene_idx(self, gene_name):
		return self.__gene_name_idx_map[gene_name]

	def test_intervals(self):
		return self.__gene_intervals
	def filter_gene_ann(self, ann):
		if "hap" in ann.chr or "_" in ann.chr:
			return True
		else:
			return False

	def load_gene_bed(self, gene_ann_bed):
		start_time = time.time()
		n = 0
		for line in open(gene_ann_bed, "r"):
			ann = GeneBed(line)
			# filter out annotation not needed, e.g. gene annotation on chr_hap
			if self.filter_gene_ann(ann):
				continue
			self.__genes.setdefault(ann.chr, []).append(ann)
			# gene id and name map, this dict will take ~3G ram for ensgene annotation, current not loaded
			#self.__gene_name_id_map.setdefault(ann.gene_name, ann.gene_id)
			n += 1
		print >> sys.stderr, n, "annotations from", gene_ann_bed, "loaded."
		# sort all ann of each chr by start pos
		for chr in self.__genes:
			self.__genes[chr] = sorted(self.__genes[chr], key=lambda d:int(d.start))

			# save start pos in self.__gene_starts
			self.__gene_starts[chr] = [int(d.start) for d in self.__genes[chr]]
		print >> sys.stderr, time.time() - start_time, "sec. elapsed."	

	# for a given gene, return its interval if it is in the annotation, return empty interval if not.
	def get_gene_interval(self, gene_name):
		if gene_name in self.__gene_intervals:
			return self.__gene_intervals[gene_name]
		else:
			print >> sys.stderr, "Warnning: gene name", gene_name, "is not in current annotation."
			return GeneIntervals([])
	# get distance between two gene intervals, return -1 if genes are rom different chrs, minus value if overlap.
	def get_gene_distance(self, gene1, gene2):
		interval1 = self.get_gene_interval(gene1)
		interval2 = self.get_gene_interval(gene2)
		if interval1 and interval2:
			chr1 = interval1[0]
			start1 = int(interval1[1])
			end1 = int(interval1[2])
			chr2 = interval2[0]
			start2 = int(interval2[1])
			end2 = int(interval2[2])
			if chr1 != chr2: # genes from different chrs
				return -1 
			else:
				return cmp_intervals(start1, end1, start2, end2)
			
			
			
	#return a list of GeneBed annotation		
	def map_pos_to_genes(self, chr, pos):
		# if pos is within 10bp window of any boundary, set close_to_boundary True
		t = 10
		matched_genes = []
		if not chr in self.__gene_starts:
			return matched_genes
		idx = bisect.bisect(self.__gene_starts[chr], pos)	
		while 0 < idx <= len(self.__gene_starts[chr]):
			#bpann is an GeneBed object
			bpann = self.__genes[chr][idx-1]
			#search within a limited region (default 1000000)
			if pos - bpann.start > self.__max_diff: 
				break

			if bpann.start <= pos <= bpann.end:
				'''
				# check if pos is on/close to current annotation's boundary, these are not gene annotations' attributes, need to reset every time
				bpann.is_on_boundary = False
				bpann.close_to_boundary = False
				if bpann.start == pos or bpann.end == pos:
					bpann.is_on_boundary = True
					bpann.close_to_boundary = True
				elif min(abs(bpann.start-pos), abs(bpann.end-pos)) < t:
					bpann.close_to_boundary = True
				'''
				matched_genes.append(bpann)
			idx -= 1
	
		return matched_genes
	def get_adjacent_exons(self, chr, pos):
		previous_exons = []
		next_exons = []	
			
		if not chr in self.__gene_starts:
			return previous_exons, next_exons
		idx = bisect.bisect(self.__gene_starts[chr], pos)
		while 0 < idx <= len(self.__gene_starts[chr]):
			
			bpann = self.__genes[chr][idx-1]
			#search within a limited region (default 1000000)
			if pos - bpann.start > self.__max_diff: 
				break
			if bpann.start <= pos <= bpann.end:
				# when breakpoint map to an intron/exon annotation
				#if bpann.type == "intron" or bpann.type == "cds": 
				if True: 
					# previous and next exons are based on coordinates, not strand. i.e. if a adjacent exon's coordinate < bp, it is a previous exon, otherwise next exon
					for i in range(idx, min(idx+100, len(self.__gene_starts[chr]))):
						adjacent_bpann = self.__genes[chr][i]
						# in next 100 annotations, try to find one exon next to current intron of the same transcript
						if bpann.transcript_id == adjacent_bpann.transcript_id and adjacent_bpann.type == "cds" and abs(adjacent_bpann.idx - bpann.idx) <= 1:

							next_exons.append(adjacent_bpann)
							#break
							#if bpann.strand == "f":
							#	next_exons.append(adjacent_bpann)
							#	break
							#else:
							#	previous_exons.append(adjacent_bpann)
							#	break
					for i in range(idx-2, max(0, idx-100), -1):
						adjacent_bpann = self.__genes[chr][i]
						# in next 100 annotations, try to find one exon next to current intron of the same transcript
						if bpann.transcript_id == adjacent_bpann.transcript_id and adjacent_bpann.type == "cds" and abs(adjacent_bpann.idx - bpann.idx) <= 1:
							previous_exons.append(adjacent_bpann)
							#break
							#if bpann.strand == "r":
							#	next_exons.append(adjacent_bpann)
							#	break
							#else:
							#	previous_exons.append(adjacent_bpann)
							#	break
			idx -= 1
		return previous_exons, next_exons	
	def get_adjacent_introns(self, chr, pos):
		previous_introns = []
		next_introns = []	
			
		if not chr in self.__gene_starts:
			return previous_introns, next_introns
		idx = bisect.bisect(self.__gene_starts[chr], pos)
		while 0 < idx <= len(self.__gene_starts[chr]):
			
			bpann = self.__genes[chr][idx-1]
			#search within a limited region (default 1000000)
			if pos - bpann.start > self.__max_diff: 
				break
			if bpann.start <= pos <= bpann.end:
				# when breakpoint map to an intron/exon annotation
				#if bpann.type == "intron" or bpann.type == "cds": 
				if True: 
					# previous and next exons are based on coordinates, not strand. i.e. if a adjacent exon's coordinate < bp, it is a previous exon, otherwise next exon
					for i in range(idx, min(idx+100, len(self.__gene_starts[chr]))):
						adjacent_bpann = self.__genes[chr][i]
						# in next 100 annotations, try to find one exon next to current intron of the same transcript
						if bpann.transcript_id == adjacent_bpann.transcript_id and adjacent_bpann.type == "intron" and abs(adjacent_bpann.idx - bpann.idx) <= 1:
							next_introns.append(adjacent_bpann)
					for i in range(idx-2, max(0, idx-100), -1):
						adjacent_bpann = self.__genes[chr][i]
						# in next 100 annotations, try to find one exon next to current intron of the same transcript
						if bpann.transcript_id == adjacent_bpann.transcript_id and adjacent_bpann.type == "intron" and abs(adjacent_bpann.idx - bpann.idx) <= 1:
							previous_introns.append(adjacent_bpann)
			idx -= 1
		return previous_introns, next_introns	
				
# Deprecated, use CffFusion instead
class GeneFusions():
	__fusions = {}
	__genes = {}
	__geneann = GeneAnnotation("")
	def __init__(self, fusion_file, gene_ann_bed):
		self.__geneann = GeneAnnotation(gene_ann_bed)
		self.load_gene_fusions(fusion_file)
		
	# fusion_file is a tsv file in cff(common fusion format) format
	'''
	1 gene_chromosome1
	2 genomic_break_pos1
	3 genomic_strand1
	4 gene_chromosome2
	5 genomic_break_pos2
	6 genomic_strand2
	7 orf
	8 read_through
	9 splitr_count
	10 span_count
	11 Sample
	12 Run
	13 Tool
	14 id
	15 probability
	'''
	#The original genes reported by each tool are not included, will map the locations to gene to uniform the gene names
	def load_gene_fusions(self, fusion_file):
		for line in open(fusion_file, "r"):
			if line.startswith(("@", "#", "Chromosome1", "gene_")):
				continue

			fusion = CffFusion(line)
			gene_info1 = self.__geneann.map_pos_to_genes(fusion.chr1, fusion.pos1)
			gene_info2 = self.__geneann.map_pos_to_genes(fusion.chr2, fusion.pos2)
			gene_name1 = sorted(list(set([g[7] for g in gene_info1])))
			gene_name2 = sorted(list(set([g[7] for g in gene_info2])))
			gene_transcripts1 = ",".join(sorted(list(set([g[3] for g in gene_info1]))))
			gene_transcripts2 = ",".join(sorted(list(set([g[3] for g in gene_info2]))))

			# all mapped genes
			fusion.gene1 = ",".join(gene_name1)
			fusion.gene2 = ",".join(gene_name2)
			
			# all transcripts of each mapped gene
			fusion.trans_id1 = gene_transcripts1
			fusion.trans_id2 = gene_transcripts2

			# use gene pair of each fusion as key, for recurrent fusions
			if len(gene_name1)>0 and  len(gene_name2)>0:
				key = fusion.gene1 + "_" + fusion.gene2
				reverse_key = fusion.gene2 + "_" + fusion.gene1
				
				if not reverse_key in self.__fusions:
					self.__fusions.setdefault(key, []).append(fusion)
				else:
					self.__fusions.setdefault(reverse_key, []).append(fusion)

				# use one gene at a time as key, for recurrent genes
				for gene in fusion.gene1.split(",") + fusion.gene2.split(","):
					self.__genes.setdefault(gene, []).append(fusion)
			else:
				# some process for breakpoints mapped to one gene or no genes
				pass
	# get a set of sample names, remove duplicate samples		
	def get_samples_and_tools(self, fusion_list):
		#print [f.sample for f in fusion_list]
		return list(set([f.sample for f in fusion_list])), list(set([f.tool for f in fusion_list]))
	# get partner genes
	def get_partner_genes(self, gene, fusion_list):
		partner_genes = list(set([f.gene1 for f in fusion_list if f.gene1!=gene] + [f.gene2 for f in fusion_list if f.gene2!=gene]))
		return partner_genes
		
	def count_sample_number(self, sample_list):
		n_tumor = 0
		n_normal = 0
		for sample in sample_list:
			if sample.endswith("N"):
				n_normal += 1
			else:
				n_tumor += 1
		return n_normal, n_tumor
	def output_recurrent_fusions(self):
		# fusion looks like CHP1,MEP_TIMP2
		for fusion in self.__fusions:
			genes = fusion.split("_")
			genes1 = set(genes[0].split(","))
			genes2 = set(genes[1].split(","))
			mapped_to_same_gene = bool(genes1 & genes2)

			gene_ids1 = []
			gene_ids2 = []
			for gene in genes1:
				gene_ids1.append(self.__geneann.get_gene_id(gene))
			for gene in genes2:
				gene_ids2.append(self.__geneann.get_gene_id(gene))
		

			sample_fusion_list= self.__fusions[fusion]
			sample_list, tool_list = self.get_samples_and_tools(sample_fusion_list)
						
			n_normal, n_tumor = self.count_sample_number(sample_list)
			print "Gene_Fusion:", fusion, "Sample_num:", len(sample_list), "Gene_ID:", ",".join(gene_ids1) + "_" + ",".join(gene_ids2),  "Same_gene:", mapped_to_same_gene, "Tools:", "|".join(tool_list), "#Normal:", n_normal, "#Tumor", n_tumor, "Samples:", "|".join(sample_list)
			for sample_fusion in sample_fusion_list:
				print "\t" + sample_fusion.print_to_string()
				
	def output_recurrent_genes(self):
		for gene in self.__genes:
			sample_fusion_list= self.__genes[gene]
			sample_list, tool_list = self.get_samples_and_tools(sample_fusion_list)
			partner_genes = self.get_partner_genes(gene, sample_fusion_list)

			n_normal, n_tumor = self.count_sample_number(sample_list)
			
			
			print "Gene:", gene, "Sample_num:", len(sample_list), "#Normal:", n_normal, "#Tumor", n_tumor,  "Samples:", "|".join(sample_list), "Partner_num:", len(partner_genes), "Partner_genes:", "|".join(partner_genes)
			for sample_fusion in sample_fusion_list:
				print "\t" + sample_fusion.print_to_string()

	def test(self):
		self.output_recurrent_fusions()
		self.output_recurrent_genes()
		
		
def cmp_fusion_breakpoints(bp1, bp2, diff):
	chr1, pos1, strand1 = bp1
	chr2, pos2, strand2 = bp2
			
	if chr1 != chr2 or strand1 != strand2:
		return False
	elif abs(pos1 - pos2) < diff:
		return True
	else:
		return False	
#compare two genomic locations, return True if the difference is samller than parameter diff		
def cmp_breakpoints(chr1, pos1, chr2, pos2, diff):
	if chr1 != chr2:
		return False
	elif abs(pos1 - pos2) < diff:
		return True
# check whether two intervals overlap, overlap if end-start > 0
def cmp_intervals(start1, end1, start2, end2):
	start = max(start1, start2)
	end = min(end1, end2)
	return end - start


def pass_filter(read):
	if read.mapping_quality == 255: # tophat2 max mapq
		return True
	else:
		return False
# return #read per basepair for all exons(cds) of transcripts in trans_id_list
def get_exon_cov(trans_id_list, gene_ann, bam):

	# get uniquely mapped read number in each cds
	for trans_id in trans_id_list:
		for transcript in  gene_ann.get_transcripts_ann(trans_id):
			read_cnt = 0
			n = 0
			#if type == "cds":
			if 1:
				for read in bam.fetch(transcript.chr, transcript.start, transcript.end):
					n += 1
					# bam.fetch may return reads not in given region (bug of pysam?)
					# also require read to pass a "filter", currently only requires mapping quality to be maximum, 50 for tophat2, 60 for bwa, 255 for star
					if read.reference_start in range(transcript.start, transcript.end) and pass_filter(read):
						read_cnt += 1
					else:
						pass
						#print read.reference_start, read.reference_start in [start, end]
			
				print transcript.tostring, read_cnt, float(read_cnt)/(transcript.end - transcript.start)
def get_fusion_exon_cov(gene_fusion, gene_ann, bam):
	# map fusion breakpoints to genes
	matched_genes1 = gene_ann.map_pos_to_genes(gene_fusion.chr1, gene_fusion.pos1)	
	matched_genes2 = gene_ann.map_pos_to_genes(gene_fusion.chr2, gene_fusion.pos2)	
	
	# get all transcripts of genes
	trans_id_list1 = list(set(g[3] for g in matched_genes1))
	trans_id_list2 = list(set(g[3] for g in matched_genes2))
	trans_id_list = trans_id_list1 + trans_id_list2
	get_exon_cov(trans_id_list, gene_ann, bam)

#build fusion reference, refs is pysam.FastaFile, seg_len is the length of extraced ref sequence
def output_fusion_fasta(fusion, refs, seg_len):
	## defuse strands have different meanings, need to convert to the following pattern before uesd
	# ++ : ->| |->
	# +- : ->| <-|
	# -+ : |<- |->
	# -- : |<- <-|
	## defuse:
	# ++ : ->| <-|
	# +- : ->| |->
	# -+ : |<- <-|
	# -- : |<- |->
	
	# pysam.FastaFile uses 0-based coordinates, cff fusion uses 1-based
	# strand used here is differnt from defuse
	if strand1 == "+": ## ->|
		win_start1 = fusion.pos1 - seg_len
		win_end1 = fusion.pos1
		seq1 = refs.fetch(fusion.chr1, win_start1, win_end1)
	else: ## |<-
		win_start1 = fusion.pos1 - 1
		win_end1 = fusion.pos1 + seg_len - 1
		seq1 = rc_seq(refs.fetch(fusion.chr1, win_start1, win_end1))

	if strand2 == "+": ## |->
		win_start2 = fusion.pos2 - 1
		win_end2 = fusion.pos2 + seg_len - 1
		seq2 = refs.fetch(fusion.chr2, win_start2, win_end2)
	else: ## <-|
		win_start2 = fusion.pos2 - seg_len
		win_end2 = fusion.pos2
		seq2 = rc_seq(refs.fetch(fusion.chr2, win_start2, win_end2))
	seq = seq1 + seq2
# build transcript sequence based on bed format annotation file
# chr1    24683494	24685032	ENST00000003583 utr3    0       r       STPG1   ENSG00000001460
	
def build_transcript_seq(gene_ann, ref):
	trans_seqs = {}
	for chr in gene_ann.__transcripts_ann:
		seq = []
		for transcript in gene_ann.__transcripts_ann[chr]: 	
			seq.append(refs.fetch(transcript.chr, transcript.start, transcript.end))		
				
		trans_seqs.setdefault(transcript, "".join(seq))
	return trans_seqs

#split gene order string to get gene names
# MYH8_cds>>MYH13_intron;RP11-401O9.3_utr5>>AK097500_utr3;CTC-297N7.11_utr5;RP11-799N11.1_utr5
def get_gene_names_from_gene_order(gene_order_string):
	# gene name list
	g_driver = []
	g_passenger = []
	tmp = gene_order_string.split(";")
	for t in tmp:
		tmp2 = t.split(">>")
		for n in range(0, len(tmp2)):
			t2 = tmp2[n]
			if n == 0:
				g = g_driver
			else:
				g = g_passenger
			tmp3 = t2.split(",")
			for t3 in tmp3:
				if "_" in t3:
					tmp4 = t3.split("_")
					g.append(tmp4[0])
	return g_driver, g_passenger
# build fusion sequence, with up/downstream genome sequence, splicing not considered by far		
def get_fusion_seq(fusion, ref, seg_len):

	chr1 = fusion.chr1
	chr2 = fusion.chr2
	bp1 = fusion.pos1
	bp2 = fusion.pos2
	strand1 = fusion.strand1
	strand2 = fusion.strand2

	refs = pysam.FastaFile(ref)
	if not chr1 in refs:
		if chr1.startswith("chr"):
			chr1 = chr1[3:]
			chr2 = chr2[3:]
		else:
			chr1 = "chr" + chr1
			chr2 = "chr" + chr2
	if not chr1 in refs:
		print >> sys.stderr, "unmatched reference."
		sys.exit(1)


	if strand1 == "+": ## ->|
		win_start1 = bp1 -seg_len
		win_end1 = bp1
		seq1 = refs.fetch(chr1, win_start1, win_end1)
	else: ## |<-
		win_start1 = bp1 - 1
		win_end1 = bp1 +seg_len - 1
		seq1 = sequtils.rc_seq(refs.fetch(chr1, win_start1, win_end1), "rc")

	if strand2 == "+": ## <-|
		win_start2 = bp2 - seg_len
		win_end2 = bp2
		seq2 = sequtils.rc_seq(refs.fetch(chr2, win_start2, win_end2), "rc")
	else: ## |->
		win_start2 = bp2 - 1
		win_end2 = bp2 + seg_len - 1
		seq2 = refs.fetch(chr2, win_start2, win_end2)
	fusion.seq1 = seq1.upper()
	fusion.seq2 = seq2.upper()
	if not fusion.seq1:
		fusion.seq1 = "NA"
	if not fusion.seq2:
		fusion.seq2 = "NA"
		
	refs.close()
	#return seq1.upper(), seq2.upper()
def format_chr(chr0, ref):	
	formatted_chr = chr0
	if not formatted_chr in ref:
		if formatted_chr.startswith("chr"):
			formatted_chr = formatted_chr[3:]
		else:
			formatted_chr = "chr" + formatted_chr
	
	return formatted_chr


# for a given gene, according to gene_order and bp get all its transcript sequences and 100bp potential fusion sequences (only for head gene)
def build_transcript_and_fusion_seq(gene_ann, gene_name, ref, bp, gene_order):
	
	# get all transcripts for head gene	
	seqs = []
	gene_intervals = gene_ann.get_gene_interval(gene_name)
	tids = gene_intervals.transcript_ids
	if not tids:
		print >> sys.stderr, "Gene not in annotation:", gene_name
		print >> sys.stderr, "Can not build ref seq"
		#sys.exit(1)
		return seqs
		

	is_head = (gene_order == "head")
	for transcript_id in tids:
		trans_ann_list = gene_ann.get_transcripts_ann(transcript_id) # all GeneBed annotations of transcript_id
		if not trans_ann_list[0].start <= bp <= trans_ann_list[-1].end: #skip the transcript that doesnt include the breakpoint
			continue
		is_fw = (trans_ann_list[0].strand == "f")
		trans_seqs = []
		fusion_seqs = []
		if is_head:
			if is_fw: # head gene on forward strand
				for ann in trans_ann_list:
					# breakpoint in current bed annotation
					if ann.start <= bp <= ann.end:
						seq = ref.fetch(format_chr(ann.chr, ref), ann.start-1, bp)
						fusion_seqs.append(seq)
						if ann.type == "intron": # assume a intron-retention
							trans_seqs.append(seq)
					if ann.type != "intron":
						seq = ref.fetch(format_chr(ann.chr, ref), ann.start-1, ann.end)
						trans_seqs.append(seq)
						if ann.end < bp:
							fusion_seqs.append(seq)
			else:
				for ann in trans_ann_list:
					# breakpoint in current bed annotation
					if ann.start <= bp <= ann.end:
						seq = ref.fetch(format_chr(ann.chr, ref), bp-1, ann.end)
						fusion_seqs.append(seq)
						if ann.type == "intron": # breakpoint in intron, cant build seq
							trans_seqs.append(seq)
						
					if ann.type != "intron": # for now do not consider intonic breakpoint
						seq = ref.fetch(format_chr(ann.chr, ref), ann.start-1, ann.end)
						trans_seqs.append(seq)
						if ann.start > bp:
							fusion_seqs.append(seq)

		else: # tail gene, do not build transcript seq
			if is_fw: # tail gene on forward strand
				for ann in trans_ann_list:
					# breakpoint in current bed annotation
					if ann.start <= bp <= ann.end:
						seq = ref.fetch(format_chr(ann.chr, ref), bp-1, ann.end)
						fusion_seqs.append(seq)
					elif ann.start > bp and ann.type != "intron": # for now do not consider intonic breakpoint
						seq = ref.fetch(format_chr(ann.chr, ref), ann.start-1, ann.end)
						fusion_seqs.append(seq)
						
			else:
				for ann in trans_ann_list:
					# breakpoint in current bed annotation
					if ann.start <= bp <= ann.end:
						seq = ref.fetch(format_chr(ann.chr, ref), ann.start-1, bp)
						fusion_seqs.append(seq)
						break
					elif ann.type != "intron": # for now do not consider intonic breakpoint
						seq = ref.fetch(format_chr(ann.chr, ref), ann.start-1, ann.end)
						fusion_seqs.append(seq)

		trans_seq = ("".join(trans_seqs)).upper()
		fusion_seq = ("".join(fusion_seqs)).upper()
			
		'''
		print transcript_id
		print trans_seq.startswith(fusion_seq)
		print len(trans_seq), len(fusion_seq)
		print "1",fusion_seq[-100:]
		fusion_bp = len(fusion_seq)
		print "3", trans_seq[max(0, fusion_bp-seq_len):fusion_bp]
		'''
		
		if fusion_seq:
			seq_len = 100
			if not is_fw:
				trans_seq = sequtils.rc_seq(trans_seq, "rc")
				fusion_seq = sequtils.rc_seq(fusion_seq, "rc")

			if is_head:
				fusion_bp = len(fusion_seq)
				fusion_seq = fusion_seq[-seq_len:]
				trans_seq1 = trans_seq[max(0, fusion_bp-seq_len):fusion_bp]
				trans_seq2 = trans_seq[fusion_bp:fusion_bp+seq_len]

			else:
				fusion_bp = len(trans_seq) - len(fusion_seq)
				fusion_seq = fusion_seq[:seq_len]
				trans_seq1 = trans_seq[max(0, fusion_bp-seq_len):fusion_bp]
				trans_seq2 = trans_seq[fusion_bp:fusion_bp+seq_len]

			seqs.append((trans_seq1, trans_seq2, fusion_seq))
	return seqs		
# for given gene annotation bed file, build junctions sequences for echo transcript
# ref is pysam.FastaFile

#chr1    24683494        24685032        ENST00000003583 utr3    0       r       STPG1   ENSG00000001460
#chr1    24685032        24685109        ENST00000003583 cds     0       r       STPG1   ENSG00000001460
#chr1    24685109        24687340        ENST00000003583 intron  0       r       STPG1   ENSG00000001460
def output_trans_seq(trans_seqs, bp_list, cur_trans_id, cur_trans_strand, junc_seq_len):
	trans_seq = "".join(trans_seqs)
	junc_id = 0
	for bp in bp_list:
		junc_id += 1
		seq1 = trans_seq[max(0, bp-junc_seq_len):bp]
		seq2 = trans_seq[bp:min(len(trans_seq), bp+junc_seq_len)]
		seq = seq1 + seq2
		if cur_trans_strand == "r":
			seq = sequtils.rc_seq(seq, "rc")
		ref_name = ">" + cur_trans_id + "_" + str(junc_id)
		print ref_name
		#print seq1, seq2
		print seq

	
def build_junction_seq_for_gene_bed(ref, gene_bed):
	cur_trans_id = ""
	cur_trans_strand = ""
	cur_trans_chr = ""
	junc_seq_len = 100
	trans_seqs = []
	bp_list = []
	for line in open(gene_bed, "r"):
		ann = GeneBed(line)
		if not format_chr(ann.chr, ref) in ref:
			continue
		if cur_trans_id != ann.transcript_id:
			#output transcript junction sequences
			output_trans_seq(trans_seqs, bp_list, cur_trans_id, cur_trans_strand, junc_seq_len)

			cur_trans_id = ann.transcript_id
			cur_trans_strand = ann.strand
			cur_trans_chr = ann.chr
			bp_list = []
			trans_seqs = []
		if ann.type != "intron":
			seq = ref.fetch(format_chr(ann.chr, ref), ann.start-1, ann.end)	
			trans_seqs.append(seq)
		else:
			bp_list.append(len("".join(trans_seqs)))
		#print ann.tostring(), ann.end - ann.start
	output_trans_seq(trans_seqs, bp_list, cur_trans_id, cur_trans_strand, junc_seq_len)
			
				
