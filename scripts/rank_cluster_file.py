#! /usr/bin/env python

import sys
import pygeneann_MetaFusion as pygeneann
import argparse

def output_filtered_list(category_list):
	for c in category_list:
		c.out()

# return fusions with exactly "tools_num" tools
def filter_tools_num(fusion_list, tools_num):
    return filter(lambda x:len(x.tools)==tools_num, fusion_list)


parser = argparse.ArgumentParser()

parser.add_argument('fusion_cluster_file', action='store', help='cluster file)')
parser.add_argument('--TP', required=False, action='store', help='TP cluster file)')
args = parser.parse_args()

category_stats = pygeneann.CategoryFusionStats(args.fusion_cluster_file)

fusion_list = category_stats.category_list
descending_fusion_list= []
# get number of tools
i=max([len(fusion.tools) for fusion in fusion_list])
#generate list sorted by tool number in descending order
while i > 0:
    fusions_num_tools = filter_tools_num(fusion_list, i)  
    descending_fusion_list += fusions_num_tools
    i -= 1
if len(sys.argv) == 4:
    category_stats_TP = pygeneann.CategoryFusionStats(args.TP)
    fusion_list_TP = category_stats_TP.category_list
    TP_FIDs=[fusion.fusion_IDs for fusion in fusion_list_TP]
    for fusion in descending_fusion_list:
        if fusion.fusion_IDs in TP_FIDs: 
            print("TRUE POSITIVE,  # of tools: " + str(len(fusion.tools)))
            fusion.out()
        else: 
            pass
            print ("FP/UNKNOWN, # of tools: " + str(len(fusion.tools)))
            fusion.out()
elif len(sys.argv) == 2:
    # output header
    #pygeneann.output_cluster_header()
    pygeneann.output_cluster_header_subset()
    for fusion in descending_fusion_list:
        fusion.out_subset()
