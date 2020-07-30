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
#print(pygeneann.GeneAnnotation.__genes)
#print(dir(gene_ann))
#chr="chr5"
#print([gene_bed.gene_name for gene_bed in gene_ann._GeneAnnotation__genes[chr]])
#print(gene_ann._GeneAnnotation__gene_starts[chr])


def choose_closest_exon(exon_list, pos, strand, ori):
    #chooses the closest exon from a list of adjacent exons either upstream (head) or downstream (tail) of the fusion breakpoint
    smallest_dist=float('inf')
    closest_exon=""
    for exon in exon_list:
      # get dist, exon.start will always be larger value than pos2
      if (strand=="+" and ori=="tail") or (strand=="-" and ori=="head"):
        dist=exon.start - pos 
      if (strand=="+" and ori=="head") or (strand=="-" and ori=="tail"):
        dist=pos - exon.start
      # check dist and re-assign exon
      if dist <= smallest_dist: closest_exon=exon
      exon_igv=str(exon.chr) + ":" + str(exon.start) + "-" + str(exon.end)
      if debug:
        print(dist)
        print(exon_igv)
    if debug: print("Closest exon: " + str(closest_exon.chr) + ":" + str(closest_exon.start) + "-" + str(closest_exon.end))
    return closest_exon
        
gene_ann = pygeneann.GeneAnnotation(ensbed)

for line in open(cff_file, "r"):
    fusion = pygeneann.CffFusion(line)
    adjacent_exons1 = gene_ann.get_closest_exon_lists(fusion.chr1, fusion.pos1)
    adjacent_exons2 = gene_ann.get_closest_exon_lists(fusion.chr2, fusion.pos2)
    #exons downstream of breakpoint on + strand
    if debug:
      print("HEAD GENE: ", fusion.t_gene1)
      print("BREAKPOINT: ", fusion.chr1, fusion.pos1, fusion.strand1) 
    #PREV HEAD 
    if fusion.strand1 == '+':
      try: 
        exon=choose_closest_exon(adjacent_exons1[0], fusion.pos1, fusion.strand1, "head")
        exon=str(exon.chr) + ":" + str(exon.start) + "-" + str(exon.end)
      except:
        exon="NA"
      fusion.splice_site1=exon
    #NEXT HEAD
    elif fusion.strand1 == '-':
      try:
        exon=choose_closest_exon(adjacent_exons1[1], fusion.pos1, fusion.strand1, "head")
        exon=str(exon.chr) + ":" + str(exon.start) + "-" + str(exon.end)
      except:
        exon="NA"
      fusion.splice_site1=exon
    else: raise Exception('Strand must be either + or -')

    if debug: 
      print("TAIL GENE: ", fusion.t_gene2)
      print("BREAKPOINT: ", fusion.chr2, fusion.pos2, fusion.strand2) 

    #PREV TAIL 
    if fusion.strand2 == '-' or debug:
      try:
        exon=choose_closest_exon(adjacent_exons2[0], fusion.pos2, fusion.strand2, "tail")
        exon=str(exon.chr) + ":" + str(exon.start) + "-" + str(exon.end)
      except:
        exon="NA"
      fusion.splice_site2=exon
    #NEXT TAIL
    elif fusion.strand2 == '+':
      try:
        exon=choose_closest_exon(adjacent_exons2[1], fusion.pos2, fusion.strand2, "tail")
        exon=str(exon.chr) + ":" + str(exon.start) + "-" + str(exon.end)
      except:
        exon="NA"
      fusion.splice_site2=exon
    else: raise Exception('Strand must be either + or -')
    print fusion.tostring() 

