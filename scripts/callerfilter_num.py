#! /usr/bin/env python

import sys
sys.path.append("/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin")
import pygeneann_reads_capture_DEV as pygeneann
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
#filtered_list = filter_tools_num(category_stats.category_list, int(args.num_tools)) 
filtered_list = filter_tools_num(category_stats.category_list, args.num_tools)
output_filtered_list(filtered_list)
#output_filtered_list(category_stats.category_list)
