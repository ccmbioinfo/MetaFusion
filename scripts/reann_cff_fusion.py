#!/usr/bin/env python
import sys
import  pygeneann_MetaFusion as pygeneann
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
#Load bed format gene annotation
gene_ann = pygeneann.GeneAnnotation(ensbed)

n = 1    
for line in open(cff_file, "r"):
                
    fusion = pygeneann.CffFusion(line)
    orig_pos = fusion.pos1 # record the pos before shift
    # ann_gene_order is an instance method of CffFusion class.  
    # Attempts to identify 5'->3' gene order
    # If gene name or gene loc doesnt exist in gene bed, gene order is more likely to be switched
    fusion.ann_gene_order(gene_ann)

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
