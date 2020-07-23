#! /usr/bin/env python

import sys
import argparse
from pygeneann import *
import os

# instantiate parser object
parser = argparse.ArgumentParser()

parser.add_argument('fusion_cluster_file', action='store', help='Fusion reannation file clustered by head/tail genes. (merged.cff.reann.dnasupp.bwafilter.<seq_len>.cluster file)')
parser.add_argument('output_dir', action='store', help='absolute path of the output directory')

args = parser.parse_args()


def generate_filtered_category_file(fusion_stats, output_file, tool=None, num=None):
    """
    Generates a category file for a specific tool

    :param fusion_stats: a ValidatedFusionStats object
    :param tool: name of tool we are selecting for
    :return: name of category file generated, same as output_file
    """
    # filter fusions
    if tool:
        filtered_list = fusion_stats.filter_tools_name(fusion_stats.category_list, tool)
    elif num:
        filtered_list = fusion_stats.filter_tools_num(fusion_stats.category_list, num)
    category_file = open(output_file, "w+")
    category_file.write("#cluster_type  gene1   gene2   max_split_cnt   max_span_cnt    sample_type     disease tools   inferred_fusion_type    gene1_on_bnd    gene1_close_to_bnd      gene2_on_bnd    gene2_close_to_bnd      dna_supp        samples chr1    breakpoint_1    chr2    breakpoint_2    gene1_candidates        gene2_candidates        gene1_strands   gene2_strands"  + "\n")
    for category_fusion in filtered_list:
        category_file.write(category_fusion.line + "\n")
    category_file.close()

    return output_file


#define cluster outdir
cluster_outdir = os.path.join(args.output_dir, "cluster_stats_files")
if not os.path.exists(cluster_outdir):
    os.mkdir(cluster_outdir)

# create CategoryFusionStats object for overall .cluster file 
total_fusion_stats = CategoryFusionStats(args.fusion_cluster_file)

# separate calls into files per caller
tools = ["defuse", "integrate", "fusionmap", "ericscript", "star_fusion"]
for tool in tools:
    #generate category file
    category_file = generate_filtered_category_file(total_fusion_stats, os.path.join(cluster_outdir, tool + ".cluster"), tool=tool)

# GENERATE CATEGORY COUNTS FOR EACH CALLER
tool_cluster_files = [os.path.join(cluster_outdir, tool + ".cluster") for tool in tools]
print(tool_cluster_files)
category_count_file = open(os.path.join(args.output_dir, "category_count_file.txt"), 'w+')
categories = ['GeneFusion', 'SameGene', 'NoDriverGene', 'ReadThrough', 'TruncatedNoncoding', 'TruncatedCoding']
category_count_file.write("\t".join(["tool"] + categories) + "\n" )

for tool in tools:
    fusion_stats = CategoryFusionStats(os.path.join(cluster_outdir, tool + ".cluster"))
    category_dict = fusion_stats.generate_category_counts()
    print(tool, "Number of fusions:", fusion_stats.num_fusions, category_dict)
    category_count_file.write("\t".join([tool] + [str(category_dict[category]) for category in categories] ) + "\n" )
category_count_file.close()
