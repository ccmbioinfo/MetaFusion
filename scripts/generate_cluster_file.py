#!/usr/bin/env python
import sys
import pandas as pd
import numpy as np
sys.path.append("/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin")
import  pygeneann_reads_capture_DEV as pygeneann
import pybedtools.bedtool as bedtools
import itertools
import sequtils
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('cff', action='store', help='CFF file')
parser.add_argument('FIDs', action='store', help='clustered FIDs file')

args = parser.parse_args()
cff=args.cff
FIDs=args.FIDs

def output_clustered_fusions(fusion_list, cluster_type):
        use_reference_fusion=1
        #use first fusion in fusion as reference fusion to account for flipped fusions for BP_Cluster
        if use_reference_fusion:
            try:
                ref_fus = [fusion for fusion in fusion_list if fusion.tool != 'defuse'][0]
            except IndexError:
                ref_fus = fusion_list[0]
            #ref_fus_bp1 = (ref_fus.chr1, ref_fus.pos1, ref_fus.strand1)
            #flip all fusions not ordered the same way as reference:
            for fusion in fusion_list:
                #TEST: remove reannotated gene name
                #fusion.reann_gene1 = fusion.t_gene1
                #fusion.reann_gene2 = fusion.t_gene2

                #fus_bp1 = (fusion.chr1, fusion.pos1, fusion.strand1)
                # check both breakpoint distance AND gene names
                #if not cmp_fusion_breakpoints(ref_fus_bp1, fus_bp1, 100000) and not (ref_fus.reann_gene1 == fusion.reann_gene1):
                if not (ref_fus.t_gene1 == fusion.t_gene1) and ( ref_fus.t_gene1 == fusion.t_gene2 or ref_fus.t_gene2 == fusion.t_gene1 ):
                    #flip fusion:
                    fusion.chr1, fusion.chr2 = fusion.chr2, fusion.chr1
                    fusion.pos1, fusion.pos2 = fusion.pos2, fusion.pos1
                    fusion.strand1, fusion.strand2 = fusion.strand2, fusion.strand1
                    fusion.t_gene1, fusion.t_gene2 = fusion.t_gene2, fusion.t_gene1
                    fusion.t_area1, fusion.t_area2 = fusion.t_area2, fusion.t_area1
                    #fusion.reann_gene1, fusion.reann_gene2 = fusion.reann_gene2, fusion.reann_gene1
                    fusion.reann_type1,fusion.reann_type2 = fusion.reann_type2,fusion.reann_type1
        #need to remove NA values from below lists, but only if there are valid gene names in them
        #gene1_list_orig = [f.reann_gene1 for f in fusion_list]
        gene1_list_orig = [f.t_gene1 for f in fusion_list]
        gene1_list = [gene for gene in gene1_list_orig if gene !="NA"]
        if len(gene1_list) < 1:
            gene1_list = gene1_list_orig
        #gene2_list_orig = [f.reann_gene2 for f in fusion_list]
        gene2_list_orig = [f.t_gene2 for f in fusion_list]
        gene2_list = [gene for gene in gene2_list_orig if gene !="NA"]
        if len(gene2_list) < 1:
            gene2_list = gene2_list_orig

        max_split_cnt = max([f.split_cnt for f in fusion_list])
        max_span_cnt = max([f.span_cnt for f in fusion_list])
        sample_list = [f.sample_name for f in fusion_list]

        # get fusion.captured_reads averages for T and N
        captured_reads_tumor = [f.captured_reads for f in fusion_list if f.sample_type == "Tumor"]
        try:
            captured_reads_tumor_mean = sum(captured_reads_tumor)/float(len(captured_reads_tumor))
        except ZeroDivisionError:
            captured_reads_tumor_mean = -1
        captured_reads_normal = [f.captured_reads for f in fusion_list if f.sample_type == "Normal"]
        try:
            captured_reads_normal_mean = sum(captured_reads_normal)/float(len(captured_reads_normal))
        except ZeroDivisionError:
            captured_reads_normal_mean = -1

        # list of Fusion IDs
        fusion_IDs = [f.fusion_id for f in fusion_list]
        disease_list = [f.disease for f in fusion_list]
        tool_list = [f.tool for f in fusion_list]
        sample_type_list = [f.sample_type for f in fusion_list]
        gene1_on_bndry = "True" in [f.gene1_on_bndry for f in fusion_list]
        gene1_close_to_bndry = "True" in [f.gene1_close_to_bndry for f in fusion_list]
        gene2_on_bndry = "True" in [f.gene2_on_bndry for f in fusion_list]
        gene2_close_to_bndry = "True" in [f.gene2_close_to_bndry for f in fusion_list]

        dna_supp_cluster_num = max([int(f.dnasupp) for f in fusion_list])

        category_list = [f.category for f in fusion_list]
        # PRIORITIZE CATEGORIES TO REMOVE MULTIPLE CATEGORIES
        if "ReadThrough" in category_list:
            category_list = ["ReadThrough"]
        elif "GeneFusion" in category_list:
            category_list = ["CodingFusion"]
        elif "TruncatedCoding" in category_list:
            category_list = ["TruncatedCoding"]
        elif "TruncatedNoncoding" in category_list:
            category_list = ["TruncatedNoncoding"]
        elif "NoDriverGene" in category_list:
            category_list = ["NoHeadGene"]

        # added 4 new fields to final .cluster file for purposes of validation 
        chr1_list = [str(f.chr1) for f in fusion_list]
        breakpoint_1_list = [str(f.pos1) for f in fusion_list]
        chr2_list = [str(f.chr2) for f in fusion_list]
        breakpoint_2_list = [str(f.pos2) for f in fusion_list]
        # print statement modified to include the 4 above new fields
        print "\t".join(map(str, [cluster_type, ",".join(list(set(gene1_list))), ",".join(list(set(gene2_list))), max_split_cnt, max_span_cnt, ",".join(list(set(sample_type_list))), ",".join(list(set(disease_list))), ",".join(list(set(tool_list))), ",".join(list(set(category_list))), gene1_on_bndry, gene1_close_to_bndry, gene2_on_bndry, gene2_close_to_bndry, dna_supp_cluster_num, ",".join(list(set(sample_list))), ",".join(list(set(chr1_list))), ",".join(list(set(breakpoint_1_list))), ",".join(list(set(chr2_list))), ",".join(list(set(breakpoint_2_list))), captured_reads_tumor_mean, captured_reads_normal_mean,",".join(list(set(fusion_IDs)))]))


# Load cff file
lines=[line for line in open(cff, "r")]
fusion=pygeneann.CffFusion(lines[0])
header=fusion.zone1_attrs + fusion.zone2_attrs + fusion.zone3_attrs + fusion.zone4_attrs
df_cff=pd.read_csv(cff, sep='\t', keep_default_na=False, index_col=False, names=header)

# load FIDs file
FID_clusters = [line for line in open(FIDs, "r")]
#FID_clusters = FID_clusters[1] 
for cluster in FID_clusters:
    if cluster.startswith('FIDs'): continue
    FID_lst = cluster.rstrip().split(",")
    #print(FID_lst)
    #print(df_cff[df_cff['fusion_id'].isin(FID_lst)][['t_gene1', 't_gene2', 'pos1', 'pos2']])
    #print(df_cff[df_cff['fusion_id'].isin(FID_lst)].iloc[0:1:].to_string(header=False))
    
    df_cluster = df_cff[df_cff['fusion_id'].isin(FID_lst)].to_csv(header=None, index=False, sep="\t").rstrip()
    #df_cff= df_cff[df_cff['fusion_id'].isin(FID_lst)].to_csv(header=None, index=False, sep="\t").rstrip()
    fusion_list = []
    for line in df_cluster.split("\n"):
        #print(line)
        #print("TEST")
        #print(pygeneann.CffFusion(line).tostring()) 
        fusion=pygeneann.CffFusion(line)
        fusion_list.append(fusion)
    output_clustered_fusions(fusion_list, "TEST")
#df = df.to_csv(header=None, index=False).strip('\n').split('\n')
#print(pygeneann.CffFusion(df_cff[df_cff['fusion_id'].isin(FID_lst)].iloc[0:1:].to_string(header=False)))



