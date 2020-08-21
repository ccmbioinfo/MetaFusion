#!/usr/bin/env python
import sys
import pandas as pd
import numpy as np
import  pygeneann_MetaFusion as pygeneann
import itertools
import sequtils
import argparse

def output_filtered_list(category_list):
    for c in category_list:
        c.out()

parser = argparse.ArgumentParser()

parser.add_argument('fusion_cluster_file', action='store', help='Fusion reannation file clustered by head/tail genes. (.reann.cluster file)')

args = parser.parse_args()

cluster=args.fusion_cluster_file
category_stats = pygeneann.CategoryFusionStats(cluster)

fusion_list = category_stats.category_list

filtered_list=[]
diff=100000
for fusion in fusion_list:
    bp1_lst=fusion.breakpoint_1
    bp1=bp1_lst[0]
    bp2_lst=fusion.breakpoint_2
    bp2=bp2_lst[0]
    #different chromosomes
    if fusion.chr1 != fusion.chr2:
        filtered_list.append(fusion)
    # breakpoints further than diff
    elif abs(bp1 - bp2) > diff or fusion.inferred_fusion_type == "CodingFusion" or fusion.inferred_fusion_type == "ReadThrough" or fusion.inferred_fusion_type == "SameGene":
        filtered_list.append(fusion)
# output header
pygeneann.output_cluster_header()
output_filtered_list(filtered_list)

