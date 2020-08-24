#!/usr/bin/env python
import sys
import pandas as pd
import numpy as np
import pygeneann_MetaFusion as pygeneann
import pybedtools.bedtool as bedtools
import itertools
import sequtils
import argparse
import re
import ast

#python add_db_hits_to_cluster.py /MetaFusion/RUNS/BT474.KPL4.MCF7.SKBR3.Aug-20-2020/final.cluster /MetaFusion/RUNS/BT474.KPL4.MCF7.SKBR3.Aug-20-2020/cluster.preds.collected.gencode_mapped.wAnnot.CANCER_FUSIONS

#PARSER
parser = argparse.ArgumentParser()
parser.add_argument('fusion_cluster_subset', action='store', help='Fusion cluster file )')
parser.add_argument('db_hit_file', action='store', help='FusionAnnotator output file subsetted for cancer DB hits')
args = parser.parse_args()
#
##INPUTS
cluster_file_subset=args.fusion_cluster_subset
db_hit_file = args.db_hit_file
#cluster="/MetaFusion/RUNS/melanoma.CML.Aug-20-2020/final.cluster"
#db_hit_file="/MetaFusion/RUNS/melanoma.CML.Aug-20-2020/cluster.preds.collected.gencode_mapped.wAnnot.CANCER_FUSIONS.head1"
#CREATE FUSION LIST
fusion_list = []
for line in open(cluster_file_subset, "r"):
    if line.startswith("#"):
        continue
    fusion = pygeneann.CategoryFusionSubset(line)
    fusion_list.append(fusion)
#category_stats = pygeneann.CategoryFusionStats(cluster)

# Read in "cluster.preds.collected.gencode_mapped.wAnnot.CANCER_FUSIONS" file
db_hit_lines=[line for line in open(db_hit_file, "r")]
db_hit_dict = {}
for line in db_hit_lines:
    #Extract databases which are hit
    line=line.split("\t")
    FIDs=line[7]
    hits=line[-1]
    start = hits.find("\"ATTS") 
    end = hits.find(']', start)
    hits = hits[start:end+1]
    hits = ast.literal_eval(hits.split(":")[1])
    #Populate dictionary
    db_hit_dict[FIDs] = hits
# output header
pygeneann.output_cluster_header_subset()
#Annotate .cluster file column "cancer_db_hits" with database hits 
for f in fusion_list: 
    try:
      FIDs=",".join(f.fusion_IDs)
      db_hits =  db_hit_dict[FIDs] 
      f.cancer_db_hits = db_hit_dict[FIDs]
      f.out_subset()
    except:
      f.out_subset()
