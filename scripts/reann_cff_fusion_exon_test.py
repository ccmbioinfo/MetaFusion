#!/usr/bin/env python
import sys
import  pygeneann_exon_test as pygeneann
import sequtils
import pysam
import argparse

debug=0


#parser = argparse.ArgumentParser()

#parser.add_argument('cff_file', action='store', help='CFF file, can be .cff or cff.reann')
#parser.add_argument('ensbed', action='store', help='Ensemble gene file')
#parser.add_argument('ref_fa', action='store', help='Reference genome file')

#args = parser.parse_args()

testing=0
if testing:
    #cff_file="/Users/mapostolides/Desktop/mugqic_tools/python-tools/fusiontools/0.1.0/bin/reann_cff_fusion_testing/cff_files_May6/MLLT10--TAF3.cff_val_RHS_sign_swapped"
    #cff_file="/Users/mapostolides/Desktop/mugqic_tools/python-tools/fusiontools/0.1.0/bin/reann_cff_fusion_testing/cff_files_May6/ANKIB1--AKAP9.cff_val_RHS_sign_swapped"
    cff_file="/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/reann_cff_fusion_testing/FOSB--AADACL2/FOSB--AADACL2.cff"
    #ensbed="/Users/mapostolides/Desktop/mugqic_tools/python-tools/fusiontools/0.1.0/bin/reann_cff_fusion_testing/ens_known_genes_MLLT10--TAF3.bed"
    #ensbed="/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/reann_cff_fusion_testing/ens_known_genes.bed"
    ensbed="/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/reann_cff_fusion_testing/ens_known_genes.MIR548H2.FOSB.AADACL2.bed"
    ref_fa="/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/reann_cff_fusion_testing/human_g1k_v37_decoy.fasta"
else:
    parser = argparse.ArgumentParser()

    parser.add_argument('cff_file', action='store', help='CFF file, can be .cff or cff.reann')
    parser.add_argument('ensbed', action='store', help='Ensemble gene file')
    parser.add_argument('ref_fa', action='store', help='Reference genome file')

    args = parser.parse_args()
    cff_file = args.cff_file
    ensbed = args.ensbed
    ref_fa=args.ref_fa
#Load bed format gene annotation, current support knowngene.bed's format, map given genomic loactions to genens, return matched gene list


gene_ann = pygeneann.GeneAnnotation(ensbed)
#print(pygeneann.GeneAnnotation.__genes)
#print(dir(gene_ann))
chr="chr5"
#print([gene_bed.gene_name for gene_bed in gene_ann._GeneAnnotation__genes[chr]])
#print(gene_ann._GeneAnnotation__gene_starts[chr])
n = 1    
for line in open(cff_file, "r"):
                
    fusion = pygeneann.CffFusion(line)
    # ann_gene_order is an instance method of CffFusion class.  
    #fusion.ann_gene_order(gene_ann)
    #adjacent_exons = gene_ann.get_adjacent_exons(fusion.chr1, fusion.pos1)
    #adjacent_exons = gene_ann.get_closest_exon(fusion.chr1, fusion.pos1)
    #pos1=65307954
    #pos1=65307956
    #pos2=17398997
    #pos2=17398999
    pos1=fusion.pos1
    pos2=fusion.pos2
    adjacent_exons1 = gene_ann.get_closest_exon(fusion.chr1, pos1)
    #print(adjacent_exons1)
    #exit(0)
    adjacent_exons2 = gene_ann.get_closest_exon(fusion.chr2, pos2)
    print("HEAD GENE: ", fusion.t_gene1)
    print("BREAKPOINT: ", fusion.chr1, pos1, fusion.strand1) 
    #PREV HEAD 
    if fusion.strand1 == '+' or debug:
      exon=adjacent_exons1[0][0]
      print("PREV EXON:", exon.gene_name, exon.idx, exon.chr, exon.start, exon.end, exon.strand)
    #NEXT HEAD
    elif fusion.strand1 == '-' or debug:
      exon=adjacent_exons1[1][0]
      print("NEXT EXON:",exon.gene_name,exon.idx,exon.chr, exon.start, exon.end, exon.strand )
    else: raise Exception('Strand must be either + or -')

    print("TAIL GENE: ", fusion.t_gene2)
    print("BREAKPOINT: ", fusion.chr2, pos2, fusion.strand2) 
    #PREV TAIL 
    if fusion.strand2 == '-' or debug:
      exon=adjacent_exons2[0][0]
      print("PREV EXON:",exon.gene_name,exon.idx, exon.chr, exon.start, exon.end, exon.strand)
    #NEXT TAIL
    elif fusion.strand2 == '+' or debug:
      exon=adjacent_exons2[1][0]
      print("NEXT EXON:",exon.gene_name, exon.idx, exon.chr, exon.start, exon.end, exon.strand )
    else: raise Exception('Strand must be either + or -')
 
    n += 1

