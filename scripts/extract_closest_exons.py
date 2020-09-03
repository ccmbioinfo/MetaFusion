#!/usr/bin/env python
import sys
import  pygeneann_MetaFusion as pygeneann
import sequtils
import pysam
import argparse

# set debug = 3 for STX16--RAE1
debug=0

if debug == 3:
  cff_file = "/MetaFusion/EXON_RUNS/outdir-Aug-27-2020/STX16--RAE1.ALL.cff.reann.NO_SEQ.55929088.arriba.1entry" 
  ensbed="/MetaFusion/test_data/cff_test/ens_known_genes.STX16--RAE1.renamed.bed"
  ref_fa="/MetaFusion/reference_files/human_g1k_v37_decoy.fasta"
#PARSER
else:
  parser = argparse.ArgumentParser()
  parser.add_argument('cff_file', action='store', help='CFF file, can be .cff or cff.reann')
  parser.add_argument('ensbed', action='store', help='Ensemble gene file')
  parser.add_argument('ref_fa', action='store', help='Reference genome file')
  args = parser.parse_args()
  cff_file = args.cff_file
  ensbed = args.ensbed
  ref_fa=args.ref_fa

def choose_closest_exon(exon_list, pos, strand, ori):
    #chooses the closest exon from a list of adjacent exons either upstream (head) or downstream (tail) of the fusion breakpoint
    smallest_dist=float('inf')
    closest_exon=""
    if debug == 1: 
        sys.stderr.write("There are " + str(len(exon_list)) + " exons in the list\n")
        sys.stderr.write(str([",".join([str(exon.chr) + ":" + str(exon.start) + "-" + str(exon.end), 
                                         exon.gene_name, exon.type]) for exon in exon_list]) + "\n")
    for exon in exon_list:
      # get dist, exon.start will always be larger value than pos2
      if (strand=="+" and ori=="tail") or (strand=="-" and ori=="head"):
        dist=exon.start - pos 
      if (strand=="+" and ori=="head") or (strand=="-" and ori=="tail"):
        dist=pos - exon.start
      # check dist and re-assign exon
      if dist <= smallest_dist: 
        closest_exon=exon
        smallest_dist = dist
      exon_igv=str(exon.chr) + ":" + str(exon.start) + "-" + str(exon.end)
    if debug == 1: 
      sys.stdout.write("Closest exon: " + str(closest_exon.chr) + ":" + str(closest_exon.start) + "-" + str(closest_exon.end) + "\n")
      if (strand=="+" and ori=="tail") or (strand=="-" and ori=="head"): print(closest_exon.start - pos)
      elif (strand=="+" and ori=="head") or (strand=="-" and ori=="tail"): print(pos - exon.start)
    return closest_exon
        
gene_ann = pygeneann.GeneAnnotation(ensbed)

for line in open(cff_file, "r"):
    fusion = pygeneann.CffFusion(line)
    adjacent_exons1 = gene_ann.get_closest_exon_lists(fusion.chr1, fusion.pos1)
    adjacent_exons2 = gene_ann.get_closest_exon_lists(fusion.chr2, fusion.pos2)
    #exons downstream of breakpoint on + strand
    if debug == 1:
      print("HEAD GENE: ", fusion.t_gene1)
      print("BREAKPOINT: ", fusion.chr1, fusion.pos1, fusion.strand1) 
    #HEAD LEFT OF BREAKPOINT (PREV) 
    if fusion.strand1 == '+':
      try: 
        exon=choose_closest_exon(adjacent_exons1[0], fusion.pos1, fusion.strand1, "head")
        exon=str(exon.chr) + ":" + str(exon.start) + "-" + str(exon.end)
      except:
        exon="NA"
      fusion.closest_exon1=exon
    #HEAD RIGHT OF BREAKPOINT (NEXT) 
    elif fusion.strand1 == '-':
      try:
        exon=choose_closest_exon(adjacent_exons1[1], fusion.pos1, fusion.strand1, "head")
        exon=str(exon.chr) + ":" + str(exon.start) + "-" + str(exon.end)
      except:
        exon="NA"
      fusion.closest_exon1=exon
    else: raise Exception('Strand must be either + or -')

    if debug == 1: 
      print("TAIL GENE: ", fusion.t_gene2)
      print("BREAKPOINT: ", fusion.chr2, fusion.pos2, fusion.strand2) 

    #TAIL LEFT OF BREAKPOINT (PREV)  
    if fusion.strand2 == '-':
      if debug: sys.stderr.write("Choosing tail gene feature on '" + fusion.strand2 + "' strand\n")
      try:
        exon=choose_closest_exon(adjacent_exons2[0], fusion.pos2, fusion.strand2, "tail")
        exon=str(exon.chr) + ":" + str(exon.start) + "-" + str(exon.end)
      except:
        exon="NA"
      fusion.closest_exon2=exon
    #TAIL RIGHT OF BREAKPOINT (NEXT)
    elif fusion.strand2 == '+':
      try:
        exon=choose_closest_exon(adjacent_exons2[1], fusion.pos2, fusion.strand2, "tail")
        exon=str(exon.chr) + ":" + str(exon.start) + "-" + str(exon.end)
      except:
        exon="NA"
      fusion.closest_exon2=exon
    else: raise Exception('Strand must be either + or -')
    print fusion.tostring() 


#SCRAP
#print(pygeneann.GeneAnnotation.__genes)
#print(dir(gene_ann))
#chr="chr5"
#print([gene_bed.gene_name for gene_bed in gene_ann._GeneAnnotation__genes[chr]])
#print(gene_ann._GeneAnnotation__gene_starts[chr])
