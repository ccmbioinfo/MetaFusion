#!/usr/bin/env python
import sys
import pysam

def rc_seq(seq, mode):
	c_seq = seq.upper()
	c_seq = c_seq.replace("A","a")
	c_seq = c_seq.replace("T","A")
	c_seq = c_seq.replace("a","T")
	c_seq = c_seq.replace("G","g")
	c_seq = c_seq.replace("C","G")
	c_seq = c_seq.replace("g","C")
	rc_seq = c_seq[::-1]
	r_seq = seq[::-1]
	if mode == "c":
		return c_seq
	elif mode == "rc":
		return rc_seq
	elif mode == "r":
		return r_seq
	else:
		print >> sys.stderr, "Unknown mode:", mode
		print >> sys.stderr, "Original sequence returned."
		
		return seq


def get_fusion_seq(fusion, gene_ann, ref, seg_len):
	
	chr1 = fusion.chr1	
	chr2 = fusion.chr2
	bp1 = fusion.pos1
	bp2 = fusion.pos2
	strand1 = fusion.strand1
	strand2 = fusion.strand2
		
	refs = pysam.FastaFile(ref)

	if strand1 == "+": ## ->|
		win_start1 = bp1 -seg_len
		win_end1 = bp1
		seq1 = refs.fetch(chr1, win_start1, win_end1)
	else: ## |<-
		win_start1 = bp1 - 1
		win_end1 = bp1 +seg_len - 1
		seq1 = rc_seq(refs.fetch(chr1, win_start1, win_end1), "rc")

	if strand2 == "+": ## <-|
		win_start2 = bp2 - seg_len
		win_end2 = bp2
		seq2 = rc_seq(refs.fetch(chr2, win_start2, win_end2), "rc")
	else: ## |->
		win_start2 = bp2 - 1
		win_end2 = bp2 + seg_len - 1
		seq2 = refs.fetch(chr2, win_start2, win_end2)
	return seq1.upper(), seq2.upper()

