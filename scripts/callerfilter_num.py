#! /usr/bin/env python
import sys
import pygeneann_MetaFusion as pygeneann
import argparse

def output_filtered_list(category_list):
	for c in category_list:
		c.out()

# return fusions with exactly "tools_num" tools
def filter_tools_num(fusion_list, tools_num):
    return filter(lambda x:len(x.tools)>=tools_num, fusion_list)

#parser
parser = argparse.ArgumentParser()
parser.add_argument('--cluster', action='store', help='cluster file)')
parser.add_argument('--num_tools', action='store', help='number of tools', type=int)
args = parser.parse_args()

#filter fusions
category_stats = pygeneann.CategoryFusionStats(args.cluster)
filtered_list = filter_tools_num(category_stats.category_list, args.num_tools)
# output header
pygeneann.output_cluster_header()
output_filtered_list(filtered_list)
