#!/usr/bin/env python
import sys
#sys.path.append("/hpf/largeprojects/ccmbio/jiangyue/DIPG_analysis_by_samples/Scripts/pygeneann/pygenefusionann")
#sys.path.append("/hpf/largeprojects/ccmbio/jiangyue/Genap_ccm/pygenefusionann/")
import pygeneann_sandbox as pygeneann
#import pygeneann
import sequtils
import pysam
import argparse

testing = 1
if not testing:
    parser = argparse.ArgumentParser()
    parser.add_argument('cff_file', action='store', help='CFF file, can be .cff or cff.reann')
    parser.add_argument('ensbed', action='store', help='Ensemble gene file')
    parser.add_argument('ref_fa', action='store', help='Reference genome file')
    args = parser.parse_args()
    cff_file = args.cff_file
    ensbed = args.ensbed
    ref_fa = args.ref_fa

if testing:
    #cff_file="/Users/mapostolides/Desktop/mugqic_tools/python-tools/fusiontools/0.1.0/bin/reann_cff_fusion_testing/cff_files_new/merged.cff_GABBR1"
    cff_file="/Users/mapostolides/Desktop/mugqic_tools/python-tools/fusiontools/0.1.0/bin/reann_cff_fusion_testing/cff_files_include_SJHGG136_D/merged.cff_include_stjude_validated_LOC727982_AC022311.1_ARHGEF33"
    #ensbed="/Users/mapostolides/Desktop/mugqic_tools/python-tools/fusiontools/0.1.0/bin/reann_cff_fusion_testing/ens_known_genes_ensID_GABBR1.bed"
    ensbed="/Users/mapostolides/Desktop/mugqic_tools/python-tools/fusiontools/0.1.0/bin/reann_cff_fusion_testing/ens_known_genes_LOC727982_AC022311.1_ARHGEF33.bed"
    ref_fa="/Users/mapostolides/Desktop/mugqic_tools/python-tools/fusiontools/0.1.0/bin/reann_cff_fusion_testing/human_g1k_v37_decoy.fasta"


#RPRD2--LAMC2
#cff_file="/Users/mapostolides/Desktop/mugqic_tools/python-tools/fusiontools/0.1.0/bin/reann_cff_fusion_testing/cff_files/merged.cff_RPRD2--LAMC2"
#cff_file="/Users/mapostolides/Desktop/mugqic_tools/python-tools/fusiontools/0.1.0/bin/reann_cff_fusion_testing/cff_files/merged.cff_RPRD2--LAMC2_nostrand"

#INSR--ATP9A
#cff_file="/Users/mapostolides/Desktop/mugqic_tools/python-tools/fusiontools/0.1.0/bin/reann_cff_fusion_testing/cff_files/merged.cff_INSR--ATP9A"

#ARID4B--NR5A2
#cff_file="/Users/mapostolides/Desktop/mugqic_tools/python-tools/fusiontools/0.1.0/bin/reann_cff_fusion_testing/cff_files/merged.cff_ARID4B--NR5A2"
#cff_file="/Users/mapostolides/Desktop/mugqic_tools/python-tools/fusiontools/0.1.0/bin/reann_cff_fusion_testing/cff_files/merged.cff_ARID4B--NR5A2_nostrand"

#AKAP9--ANKIB1
#cff_file="/Users/mapostolides/Desktop/mugqic_tools/python-tools/fusiontools/0.1.0/bin/reann_cff_fusion_testing/cff_files/merged.cff_AKAP9--ANKIB1"

#BCL3--CTB-171A8.1
#cff_file="/Users/mapostolides/Desktop/mugqic_tools/python-tools/fusiontools/0.1.0/bin/reann_cff_fusion_testing/cff_files/merged.cff_BCL3--CTB-171A8.1"
#cff_file="/Users/mapostolides/Desktop/mugqic_tools/python-tools/fusiontools/0.1.0/bin/reann_cff_fusion_testing/cff_files/merged.cff_BCL3--CTB-171A8.1_nostrand"

#RPTOR--PRAME
#cff_file="/Users/mapostolides/Desktop/mugqic_tools/python-tools/fusiontools/0.1.0/bin/reann_cff_fusion_testing/cff_files/merged.cff_RPTOR--PRAME"

#cff_file="/Users/mapostolides/Desktop/mugqic_tools/python-tools/fusiontools/0.1.0/bin/reann_cff_fusion_testing/cff_files/merged.cff"
#ensbed="/Users/mapostolides/Desktop/mugqic_tools/python-tools/fusiontools/0.1.0/bin/reann_cff_fusion_testing/ens_known_genes_ensID.bed"
#ensbed="/Users/mapostolides/Desktop/mugqic_tools/python-tools/fusiontools/0.1.0/bin/reann_cff_fusion_testing/ens_known_genes_ensID_ARID4B--NR5A2.bed"
#ensbed="/Users/mapostolides/Desktop/mugqic_tools/python-tools/fusiontools/0.1.0/bin/reann_cff_fusion_testing/ens_known_genes_ensID_BCL3--CTB-171A8.1.bed"
#ensbed="/Users/mapostolides/Desktop/mugqic_tools/python-tools/fusiontools/0.1.0/bin/reann_cff_fusion_testing/ens_known_genes_ensID_RPRD2_LAMC2.bed"
#ensbed="/Users/mapostolides/Desktop/mugqic_tools/python-tools/fusiontools/0.1.0/bin/reann_cff_fusion_testing/ens_known_genes_ensID_AKAP9--ANKIB1.bed"
#ensbed="/Users/mapostolides/Desktop/mugqic_tools/python-tools/fusiontools/0.1.0/bin/reann_cff_fusion_testing/ens_known_genes_ensID_RPTOR--PRAME.bed"
#ref_fa="/Users/mapostolides/Desktop/mugqic_tools/python-tools/fusiontools/0.1.0/bin/reann_cff_fusion_testing/human_g1k_v37_decoy.fasta"

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

    # the "check_codon" function somehow shifts breakpoints. I've decided to omit this for now
    fusion.check_codon(gene_ann, ref_fa)
    if orig_pos != fusion.pos1: # the breakpoints have been shifted
        fusion.ann_gene_order(gene_ann) # remap genes with shifted breakpoints
    pygeneann.get_fusion_seq(fusion, ref_fa, 100)
    print fusion.tostring()
    n += 1




# SCRAP
#gene_ann2 = pygeneann.GeneAnnotation(args.ensbed)
#testing = 1
#if testing == 1:
#    pickle_file = "/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/reann_cff_fusion_testing/PICKLE_GeneAnnotation_head"
#    try:
#        fh = open(pickle_file, "rb")
#        test_pickle = pickle.load(fh)        
#        print >> sys.stderr, test_pickle
#        print >> sys.stderr, "pickled test object loaded"
#        fh.close()
#        print >> sys.stderr, "fh closed (try)"
#    except (OSError, IOError) as e:
#        test_pickle = pygeneann.GeneAnnotation(args.ensbed)
#        fh = open(pickle_file, "wb")
#        pickle.dump(test_pickle, fh)
#        print >> sys.stderr, test_pickle
#        print >> sys.stderr, "pickled test object dumped"
#        fh.close()
#        print >> sys.stderr, "fh closed (except)"
#
#
##GeneAnnotation object contains information for genomic features of each gene
#elif testing == 2:
#    pickle_file = "/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/reann_cff_fusion_testing/PICKLE_GeneAnnotation"
#    try:
#        fh = open(pickle_file, "rb")
#        gene_ann = pickle.load(fh)
#        print >> sys.stderr, "pickled GeneAnnotation object loaded"
#        fh.close()
#        print >> sys.stderr, "fh closed (try)"
#    except (OSError, IOError) as e:
#        gene_ann = pygeneann.GeneAnnotation(args.ensbed)
#        fh = open(pickle_file, "wb")
#        pickle.dump(gene_ann, fh)
#        print >> sys.stderr, "pickled GeneAnnotation object dumped" 
#        fh.close()
#        print >> sys.stderr, "fh closed (except)"
#else:
#    gene_ann = pygeneann.GeneAnnotation(args.ensbed)
#

