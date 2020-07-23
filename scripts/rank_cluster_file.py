#! /usr/bin/env python

import sys
#sys.path.append("/hpf/largeprojects/ccmbio/jiangyue/DIPG_analysis_by_samples/Scripts/pygeneann/pygenefusionann")
#sys.path.append("/hpf/largeprojects/ccmbio/jiangyue/Genap_ccm/pygenefusionann/")
sys.path.append("/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin")
import pygeneann_reads_capture_DEV as pygeneann
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
# len(sys.argv) will be 4 with TP file, and 2 without
#print(len(sys.argv))

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
#output_filtered_list(descending_fusion_list)
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
    for fusion in descending_fusion_list:
        #print("# of tools: " + str(len(fusion.tools)))
        fusion.out()
#COMMANDS/SCRAP
#python rank_cluster_file.py  /hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/testing_pipeline/BT474.KPL4.MCF7.SKBR3.benchmark.June-29-2020/merged.cff.renamed.reann.cluster.RT_filter.callerfilter2.blck_filter.ANC_filter  /hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/testing_pipeline/BT474.KPL4.MCF7.SKBR3.benchmark.June-29-2020/merged.cff.renamed.reann.cluster.RT_filter.callerfilter2.blck_filter.ANC_filter.TP 
#/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/testing_pipeline/UHRR.benchmark.June-19-2020/merged.cff.UHR.renamed.reann.cluster.RT_filter.callerfilter2.blck_filter.ANC_filter
