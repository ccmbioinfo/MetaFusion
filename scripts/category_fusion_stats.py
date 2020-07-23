#! /usr/bin/env python

import sys
#sys.path.append("/hpf/largeprojects/ccmbio/jiangyue/DIPG_analysis_by_samples/Scripts/pygeneann/pygenefusionann")
#sys.path.append("/hpf/largeprojects/ccmbio/jiangyue/Genap_ccm/pygenefusionann/")
import pygeneann_OLD_star_fusion_defuse_style as pygeneann
import argparse
import itertools
import subprocess
# instantiate parser object
parser = argparse.ArgumentParser()

parser.add_argument('fusion_cluster_file', action='store', help='Fusion reannation file clustered by head/tail genes. (.reann.cluster file)')
parser.add_argument('sample_name_file', action='store', help='A file containing all sample names')

args = parser.parse_args()

class Sampleinfo():
    def __init__(self, line):
        self.sample_name, self.disease, self.sample_type = line.split()

def output_filtered_list(category_list):
    for c in category_list:
        #print c.gene1_on_bnd
        #print c.gene1_close_to_bnd
        #print c.gene2_on_bnd
        #print c.gene2_close_to_bnd
        c.out()

def output_cnt_table(fusion_list, group_tag="NA"):
    sample_type_list = ["Tumor", "Normal", "Tumor,Normal"]
    category_list = ["ReadThrough", "GeneFusion", "TruncatedCoding", "TruncatedNoncoding", "NoDriverGene", "SameGene"]
    print "Group:", group_tag
    print "Category\t" + "\t".join(sample_type_list)
    for category in category_list:
        print category ,
        for sample_type in sample_type_list:
            filtered_list = category_stats.filter_sample_type(fusion_list, sample_type)
            filtered_list = category_stats.filter_inferred_type(filtered_list, category)
            print len(filtered_list),
        print

def get_sample_type(sample):
    if sample.startswith("DIPG"):
        if sample.endswith("T"):
            sample_type = "Tumor"
        elif sample.endswith("N"):
            sample_type = "Normal"
    elif sample.startswith("SJH"):
            sample_type = "Tumor"
    elif sample.startswith("GTEX"):
        sample_type = "Normal"
    else:
        print >> sys.stderr, "Unkonwn sample:", sample
    return sample_type

def output_sample_fusion_cnt(fusion_list, sampleinfo_dict, group_tag):
    """
    :param fusion_list: list of CategoryFusion objects
    :param sampleinfo_dict: dictionary of (sample_name : Sampleinfo object)
    :param group_tag: String representing name of sample group
    """
    sample_list = []
    for fu in fusion_list:
        sample_list.extend(fu.samples)
    out_list = []
    num = 0
    for sample in sampleinfo_dict:
        # e.g. LIS_S8_L001	TP	Total	19
        out_list.append([sample, sampleinfo_dict[sample].sample_type, group_tag, sample_list.count(sample)])
    for li in sorted(out_list, key=lambda x:(x[1],x[0])):
        print "\t".join(map(str, li))

def init_sampleinfo_objects(sample_name_file):
    """
    Initializes a dictionary of (sample_name : Sampleinfo object) 
    """
    sampleinfo_dict = {}
    for line in open(sample_name_file, 'r'):
        sampleinfo = Sampleinfo(line)
        sampleinfo_dict.setdefault(sampleinfo.sample_name, sampleinfo)
    return sampleinfo_dict

def get_tools(category_list):
    """
    Returns a list of all fusion caller tools that detected fusions
    in category_list

    :param category_list: a list of CategoryFusions objects
    """
    tools = []
    for fusion in category_stats.category_list:
        for tool in fusion.tools:
            if tool not in tools:
                tools.append(tool)
    return tools

def get_fusion_types(category_list):
    """
    Returns a list of all categories under which fusions in category_list fall

    :param category_list: a list of CategoryFusions objects
    """
    fusion_types = []
    for fusion in category_stats.category_list:
        if fusion.inferred_fusion_type not in fusion_types:
            fusion_types.append(fusion.inferred_fusion_type)
    return fusion_types

def get_num_fusions_per_tool(category_stats):
    tools = get_tools(category_stats.category_list)
    fusion_counts = []
    print "tool" + "\t" + "num_detected_fusions"
    for tool in tools:
        # get filtered_list detected by tool
        filtered_list = category_stats.filter_tools_name(category_stats.category_list, tool)
        #fusion_counts.append( (tool, len(filtered_list) ) )    
        print tool + "\t" + str(len(filtered_list))
    print

def get_num_fusions_per_category(category_stats):
    categories = get_fusion_types(category_stats.category_list)
    fusion_counts = []
    print "category" + "\t" + "num_detected_fusions"
    for category in categories:
        filtered_list = category_stats.filter_inferred_type(category_stats.category_list, category)
        print category + "\t" + str(len(filtered_list)) 
    print

# init Sampleinfo object using sample_name_files
sampleinfo_dict = init_sampleinfo_objects(args.sample_name_file)

category_stats = pygeneann.CategoryFusionStats(args.fusion_cluster_file)


#get count of fusions detected by each tool
get_num_fusions_per_tool(category_stats)

# get count of fusions falling under each category
get_num_fusions_per_category(category_stats)

# output number of fusions detected in each sample
#output_sample_fusion_cnt(category_stats.category_list, sampleinfo_dict, group)

fusion_list = category_stats.category_list

# Apply Covfilter5
covfilter3=True
if covfilter3:
    fusion_list = category_stats.filter_split_cnt(fusion_list, 5)
    # need to NOT filter out FUSIONMAP calls based on span_cnt, calls are all -1, since fusionmap doesn't give this info
    fusion_list = [fusion for fusion in fusion_list if (fusion.max_span_cnt >= 5) or ('fusionmap' in fusion.tools and len(fusion.tools)==1)]


# output count table for total
print >> sys.stderr, "Total input category number:", len(fusion_list)
group = "Total"
output_cnt_table(fusion_list, group)

# enumerate all combinations of tools
print("ENUMERATING TOOLS")
combinations = []
for r in range(2,6):
    combinations += list(itertools.combinations(["ericscript", "fusionmap", "defuse", "integrate", "star_fusion"], r))
combinations = [list(combination) for combination in combinations]
print(combinations)
print(len(combinations))

#filtering based on fusion being detected by a specific list of tools, and only those tools
print "FILTERING BY TOOL COMBINATIONS AND WRITING FILE"
filename = "tool_combination_counts.tsv"
file1 = open(filename, "w+")
tool_combinations = [] # [0]
corresponding_counts = [] # [1]
for combination in combinations:
    # filter out fusions called by a given combination
    filtered_list = category_stats.filter_tools_name_multi(fusion_list, combination)
    to_write = str(sorted(combination)) + "\t" + str(len(filtered_list)) + "\n"
    tool_combinations.append(sorted(combination))
    corresponding_counts.append(len(filtered_list))
    file1.write(to_write)
file1.close()
print(tool_combinations)
#print(type(tuple(tool_combinations)))
print(corresponding_counts)

# make dictionary of {tool_combinations:corresponding_fusions}
tool_count_dict={}
for i in range(len(tool_combinations)):
    tool_count_dict[",".join(tool_combinations[i])]=corresponding_counts[i]
#tool_count_dict = dict(itertools.izip(tuple(tool_combinations),tuple(corresponding_counts)))
print(tool_count_dict)
#print "RUNNING R VENN DIAGRAM SCRIPT"
#subprocess.call ("pwd", shell=True)
#subprocess.call ("cat tool_combination_counts.tsv", shell=True)
#subprocess.call ("./make_venn_diagram.r", shell=True)
#print "R SUBPROCESS COMPLETE"



#SCRAP BELOW

#filter_sample_list = sampleinfo_dict.keys()
#filtered_list = category_stats.filter_tools_name(filtered_list, "stjude_method_Valid")
#filtered_list = category_stats.filter_samples(filtered_list, filter_sample_list)
#filtered_list = category_stats.filter_inferred_type(filtered_list, "ReadThrough")
#filtered_list = category_stats.filter_inferred_type(filtered_list, "GeneFusion")
#filtered_list = category_stats.filter_tools_name(filtered_list, "ericscript")
#output_sample_fusion_cnt(filtered_list, filter_sample_list, group)
#output_filtered_list(filtered_list)

#filtered_list = category_stats.filter_tools_num(filtered_list, 2)
#group = "2_tools"
#output_cnt_table(filtered_list, group)

#filtered_list = category_stats.filter_split_cnt(filtered_list, 5)
#filtered_list = category_stats.filter_span_cnt(filtered_list, 5)
#group = "split_span_5"
#output_cnt_table(filtered_list, group)
#output_sample_fusion_cnt(filtered_list, sampleinfo_dict, group)
#output_filtered_list(filtered_list)

"""
filtered_list = category_stats.filter_tools_num(filtered_list, 3)
group = "3_tools"
output_cnt_table(filtered_list, group)
#output_sample_fusion_cnt(filtered_list, filter_sample_list, group)

filtered_list = category_stats.filter_close_to_bnd(filtered_list)
group = "on_bnd"
output_cnt_table(filtered_list, group)
#output_sample_fusion_cnt(filtered_list, filter_sample_list, group)

filtered_list = category_stats.filter_dna_supp(filtered_list)
group = "dna_supp"
output_cnt_table(filtered_list, group)
#output_sample_fusion_cnt(filtered_list, filter_sample_list, group)
filtered_list = category_stats.filter_split_cnt(filtered_list, 20)
filtered_list = category_stats.filter_span_cnt(filtered_list, 20)
group = "split_span_20"
output_cnt_table(filtered_list, group)
#output_sample_fusion_cnt(filtered_list, args.sample_name_file, group)
"""

"""
Sample filters
"""
#filtered_list = category_stats.filter_recurrent(filter_list, 5)
#print filtered_list[0].disease
#filtered_list = category_stats.filter_disease(filtered_list, "DIPG")
#filtered_list = category_stats.filter_sample_type(filtered_list, "Tumor")
#filtered_list = category_stats.filter_sample_number(filtered_list, 3, "DIPG")
#filtered_list = category_stats.filter_split_cnt(filtered_list, 10)
#filtered_list = category_stats.filter_span_cnt(filtered_list, 10)
#filtered_list = category_stats.filter_tools_name(filtered_list, "defuse")
#filtered_list = category_stats.filter_tools_num(filtered_list, 2)
#filtered_list = category_stats.filter_dna_supp(filtered_list)
#filtered_list = category_stats.filter_inferred_type(filtered_list, "TruncatedCoding")
#filtered_list = category_stats.filter_close_to_bnd(filtered_list)
#filtered_list = category_stats.filter_recurrent(5)

#output_cnt_table(filtered_list)
#output_sample_fusion_cnt(filtered_list, args.sample_name_file)
#filtered_list = category_stats.filter_recurrent(filtered_list, 2)

#output_filtered_list(filtered_list)



