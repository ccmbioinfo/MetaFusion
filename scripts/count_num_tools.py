#!/usr/bin/env python
import sys
import os
import pandas as pd
import numpy as np
sys.path.append("/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin")
import  pygeneann_reads_capture_DEV as pygeneann
import itertools
import sequtils
import argparse


def output_filtered_list(category_list):
    for c in category_list:
        c.out()

parser = argparse.ArgumentParser()

parser.add_argument('fusion_cluster_file', action='store', help='Fusion reannation file clustered by head/tail genes. (.reann.cluster file)')
parser.add_argument('outdir', action='store', help='outdir to write files to')

args = parser.parse_args()

cluster=args.fusion_cluster_file
outdir=args.outdir
category_stats = pygeneann.CategoryFusionStats(cluster)

fusion_list = category_stats.category_list

for i in range(2,8):
  filtered_list = filter(lambda x:len(x.tools)==i, fusion_list) 
  #print(i, len(filtered_list))
  #if i == 5:
  f= open(os.path.join(outdir, cluster + "." + str(i) + "CALLERS"),"w+")
  sys.stdout = f
  output_filtered_list(filtered_list)
  f.close()

