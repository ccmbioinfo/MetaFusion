#!/usr/bin/env python
import sys
#sys.path.append("/hpf/largeprojects/ccmbio/jiangyue/DIPG_analysis_by_samples/Scripts/pygeneann/pygenefusionann")
#sys.path.append("/hpf/largeprojects/ccmbio/jiangyue/Genap_ccm/pygenefusionann/")
#import pygeneann_OLD_star_fusion_defuse_style as pygeneann
import  pygeneann_reads_capture_DEV as pygeneann
#import pygeneann_STEPH as pygeneann
import sequtils
import pysam
import argparse



parser = argparse.ArgumentParser()

parser.add_argument('--cff', action='store', help='CFF file, can be .cff or cff.reann')
parser.add_argument('--gene_bed', action='store', help='Ensemble gene file')
parser.add_argument('--ref_fa', required=False, action='store', help='Reference genome file')

args = parser.parse_args()
cff_file = args.cff
ensbed = args.gene_bed
# Assign reference fasta if provided by user
if args.ref_fa is not None:
  ref_fa=args.ref_fa
#print(args)
#print(args.ref_fa is not None )
#exit(0)
#Load bed format gene annotation, current support knowngene.bed's format, map given genomic loactions to genens, return matched gene list
gene_ann = pygeneann.GeneAnnotation(ensbed)

n = 1    
for line in open(cff_file, "r"):
                
    fusion = pygeneann.CffFusion(line)
    orig_pos = fusion.pos1 # record the pos before shift
    # ann_gene_order is an instance method of CffFusion class.  
    fusion.ann_gene_order(gene_ann)

    # TEST...
#    matched_genes1 = gene_ann.map_pos_to_genes(fusion.chr1, fusion.pos1)
#    print >> sys.stderr, matched_genes1
#    print >> sys.stderr, gene_ann == gene_ann2
    # ... TEST


    #annotate fusion id and seq
    fusion.fusion_id = "F" + (str(n)).zfill(8)

    #fusion.check_codon(gene_ann, ref_fa)
    #if orig_pos != fusion.pos1: # the breakpoints have been shifted
    #    fusion.ann_gene_order(gene_ann) # remap genes with shifted breakpoints

    # Get fusion seq only if specified by user
    if args.ref_fa is not None: 
      pygeneann.get_fusion_seq(fusion, ref_fa, 100)
    print fusion.tostring()
    n += 1

