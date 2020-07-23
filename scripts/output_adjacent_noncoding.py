#!/usr/bin/env python
import sys
import pandas as pd
import numpy as np
sys.path.append("/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin")
import  pygeneann_reads_capture_DEV as pygeneann
import itertools
import sequtils
import argparse

#cluster="/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/testing_pipeline/NEG_CONTROL_BEERS.benchmark.June-5-2020/merged.cff.renamed.reann.cluster.gene_names_cleaned.RT_filter.blck_filter.callerfilter2.FP"

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
    # breakpoints closer than diff
    if fusion.chr1 == fusion.chr2  and abs(bp1 - bp2) < diff and fusion.inferred_fusion_type != "CodingFusion" and fusion.inferred_fusion_type != "ReadThrough" and fusion.inferred_fusion_type != "SameGene":
        filtered_list.append(fusion)

output_filtered_list(filtered_list)

