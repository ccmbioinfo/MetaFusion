#! /usr/bin/env python

import sys
import argparse
from pygeneann import *
import os
import sqlite3
#from validation_classes import *
import itertools

#import numpy as np
#import pandas as pd

# instantiate parser object
parser = argparse.ArgumentParser()

parser.add_argument('validated_fusions_file', action='store', help='A file containing a set of validated fusions')
parser.add_argument('fusion_cluster_file', action='store', help='Fusion reannation file clustered by head/tail genes. (merged.cff.reann.dnasupp.bwafilter.<seq_len>.cluster file)')
#parser.add_argument('sample_name', action='store', help='Name of the sample being validated')
parser.add_argument('output_dir', action='store', help='absolute path of the output directory')

args = parser.parse_args()

class ValidatedFusionStats(CategoryFusionStats):
    """
    Takes as input a fusion_cluster_file/category_file, and validated_fusions_file.
    Generates statistics about false positives, false negatives, and true positives with respect to
    the validated fusions present in the validated_fusions_file.
    """

    def __init__(self, cluster_outdir, validated_fusions_file, fusion_cluster_file, name):
        # cluster file is loaded into self.category_list in CategoryFusionStats class definition
        CategoryFusionStats.__init__(self, fusion_cluster_file)
        # instantiate instance variables
        self.name = name
        self.cluster_files_dir = cluster_outdir
        # list of true positive CategoryFusion objects
        self.category_true_positives = []
        # list of ValidatedFusion objects
        self.validated_fusions = []
        self.true_positive_validated_fusions = []
        self.false_negative_validated_fusions = []
        # lists of gene pair tuples       
        self.false_positives = []
        self.true_positive_gene_pairs = []
        # counts
        self.num_true_positives = 0
        self.num_unvalidated_fusions = 0 
        self.num_false_negatives = 0 
        #make cluster_outdir directory if it doesn't already exist
        if not os.path.exists(self.cluster_files_dir):
            os.mkdir(self.cluster_files_dir)
        self.__load_validated_fusions_file(validated_fusions_file)
        # load category file, using altrenate __load

    def __load_validated_fusions_file(self, validated_fusions_file):
        validated_fusions_file = open(validated_fusions_file, "r")
        for line in validated_fusions_file:
            if line.startswith('#'): continue
            fusion = ValidatedFusion(line)
            self.validated_fusions.append(fusion)        
        validated_fusions_file.close()
        #print "number of validated fusions: ", len(self.validated_fusions) 

    def find_true_positives_using_candidate_gene_list(self):
        """
        Intersection of validated_fusions and output_fusions using candidate gene list obtained from genomic feature file
        These are all the validated fusions detected by the pipeline
        """
        true_positives_cluster_file_name = os.path.join(self.cluster_files_dir, self.name + ".true_positives.cluster" )
        true_positives_cluster_file = open(true_positives_cluster_file_name, "w+")
        true_positives_cluster_file.write("#cluster_type	gene1	gene2	max_split_cnt	max_span_cnt	sample_type	disease	tools	inferred_fusion_type	gene1_on_bnd	gene1_close_to_bnd	gene2_on_bnd	gene2_close_to_bnd	dna_supp	samples	chr1	breakpoint_1	chr2	breakpoint_2	gene1_candidates	gene2_candidates	gene1_strands	gene2_strands"  + "\n")
        # unvalidated fusions
        self.unvalidated_fusions = self.category_list
        for fusion in self.validated_fusions:
            for category_fusion in self.category_list:
                category_gene_pair = (category_fusion.gene1, category_fusion.gene2)
                # check if fusion candidate gene list is the samae as validated fusion candidate gene list, and filter out duplicate calls
                fusions_match = (fusion.gene1_candidates == category_fusion.gene1_candidates and fusion.gene2_candidates == category_fusion.gene2_candidates)
                # account for defuse inverted gene order calls
                fusions_match_reversed_order = (fusion.gene1_candidates == category_fusion.gene2_candidates and fusion.gene2_candidates == category_fusion.gene1_candidates)
                #samples_match = (fusion.sample_name in category_fusion.samples)
                samples_match = True 
                if ( (fusions_match or fusions_match_reversed_order) and samples_match): 
                    # flip gene order
                    if fusions_match_reversed_order:
                        category_gene_pair = (category_gene_pair[1], category_gene_pair[0])
                    if category_gene_pair not in self.true_positive_gene_pairs: 
                        # append fusion to true positive list
                        self.category_true_positives.append(category_fusion)
                        self.true_positive_gene_pairs.append(category_gene_pair)
                        true_positives_cluster_file.write(category_fusion.line + "\n") 
                    # remove  true positives, including duplicates,  to generate unvalidated_fusions list
                    self.unvalidated_fusions.remove(category_fusion)

                    if fusion not in self.true_positive_validated_fusions:
                        self.true_positive_validated_fusions.append(fusion)
        true_positives_cluster_file.close()
        self.num_true_positives = len(self.true_positive_gene_pairs)
        return true_positives_cluster_file_name

    def find_true_positives_using_breakpoints(self):
        """
        Intersection of validated_fusions and output_fusions using breakpoints 
        These are all the validated fusions detected by the pipeline
        """
        true_positives_cluster_file_name = os.path.join(self.cluster_files_dir, self.name + ".true_positives.cluster" )
        true_positives_cluster_file = open(true_positives_cluster_file_name, "w+")
        true_positives_cluster_file.write("#cluster_type	gene1	gene2	max_split_cnt	max_span_cnt	sample_type	disease	tools	inferred_fusion_type	gene1_on_bnd	gene1_close_to_bnd	gene2_on_bnd	gene2_close_to_bnd	dna_supp	samples	chr1	breakpoint_1	chr2	breakpoint_2	gene1_candidates	gene2_candidates	gene1_strands	gene2_strands"  + "\n")
        # unvalidated fusions
        self.unvalidated_fusions = self.category_list
        # define "slop" on breakpoints
        slop = 10
        for fusion in self.validated_fusions:
            for category_fusion in self.category_list:
                category_gene_pair = (category_fusion.gene1, category_fusion.gene2)
                # apply slop to multiple breakpoints for a given gene
                bp1_list = list(itertools.chain(*[range(bp - slop, bp + slop) for bp in category_fusion.breakpoint_1]))
                bp2_list = list(itertools.chain(*[range(bp - slop, bp + slop) for bp in category_fusion.breakpoint_2]))
                # check if fusion breakpoints are within genomic loci and filter out duplicate calls
                #if ( (fusion.start1 <= min(category_fusion.breakpoint_1) and min(category_fusion.breakpoint_1) <= fusion.end1) and 
                #    ( fusion.start2 <= max(category_fusion.breakpoint_2) and max(category_fusion.breakpoint_2) <= fusion.end2) ): 
                #if ( ( fusion.start1 in range(min(category_fusion.breakpoint_1) - slop, max(category_fusion.breakpoint_1) + slop ) ) and 
                #    ( fusion.start2 in range(min(category_fusion.breakpoint_2) - slop, max(category_fusion.breakpoint_2) + slop) ) ): 
                if fusion.start1 in bp1_list and fusion.start2 in bp2_list: 
                    if category_gene_pair not in self.true_positive_gene_pairs: 
                        # append fusion to true positive list
                        self.category_true_positives.append(category_fusion)
                        self.true_positive_gene_pairs.append(category_gene_pair)
                        true_positives_cluster_file.write(category_fusion.line + "\n") 
                    # remove  true positives, including duplicates,  to generate unvalidated_fusions list
                    self.unvalidated_fusions.remove(category_fusion)

                    if fusion not in self.true_positive_validated_fusions:
                        self.true_positive_validated_fusions.append(fusion)
        true_positives_cluster_file.close()
        self.num_true_positives = len(self.true_positive_gene_pairs)
        return true_positives_cluster_file_name

    def find_unvalidated_fusions(self):
        # A list of CategoryFusion objects that are unvalidated detected by pipeline
        # False positives and possible true fusions that have not been validated
        unvalidated_fusions_file_name= os.path.join( self.cluster_files_dir, self.name + ".unvalidated_fusions.cluster" )
        unvalidated_fusions_file = open(unvalidated_fusions_file_name, "w+")
        unvalidated_fusions_file.write("#cluster_type	gene1	gene2	max_split_cnt	max_span_cnt	sample_type	disease	tools	inferred_fusion_type	gene1_on_bnd	gene1_close_to_bnd	gene2_on_bnd	gene2_close_to_bnd	dna_supp	samples	chr1	breakpoint_1	chr2	breakpoint_2	gene1_candidates	gene2_candidates	gene1_strands	gene2_strands"  + "\n")
        # write entries to file
        for fusion in self.unvalidated_fusions: 
                unvalidated_fusions_file.write(fusion.line + "\n")
        unvalidated_fusions_file.close()
        self.num_unvalidated_fusions = len(self.unvalidated_fusions) 

    def find_false_negatives(self):
        """
        Set difference of validated_fusions - true_positives gives false negatives 
        Finds gene fusions which are in the provided validated_fusions file, but are not
        detected by the pipeline.
        """
        self.false_negative_validated_fusions = [validated_fusion for validated_fusion in self.validated_fusions if validated_fusion not in self.true_positive_validated_fusions]  
        false_negatives_file_name = os.path.join( self.cluster_files_dir, self.name + ".false_negatives.txt" ) 
        false_negatives_file = open(false_negatives_file_name, "w+")
        false_negatives_file.write("gene1	gene2	strand_1	strand_2	chr_1	chr_2	bp_1_start	bp_1_end	bp_2_start	bp_2_end" + "\n") 
        for fusion in self.false_negative_validated_fusions:
            false_negatives_file.write(fusion.line)
        false_negatives_file.close()
        self.num_false_negatives = len(self.false_negative_validated_fusions)
 
    def calculate_sensitivity(self):
        self.sensitivity = float(self.num_true_positives)/float(len(self.validated_fusions))

    def calculate_precision(self):
        try:
            self.precision = float(self.num_true_positives)/float(self.num_true_positives + self.num_unvalidated_fusions)
        except ZeroDivisionError:
            self.precision = 0
    def compare_validated_and_output_fusions(self):
        """
        Compares the fusion pairs in the validated fusion file with the detected fusions in the outputted fusion cluster file
        Outputs results/statistics to an output file
        :param self.validated_fusions: a list of tuples containing fusion gene pairs we would expect to be detected 
                                  by the pipeline with the given input sequence data
        :param self.category_list: file merged.cff.reann.dnasupp.bwafilter.30.cluster, outputted by the pipeline
        """
        true_positives_cluster_file_name = self.find_true_positives_using_candidate_gene_list()
        #true_positives_cluster_file_name = self.find_true_positives_using_breakpoints()
        
        self.find_unvalidated_fusions()
    
        self.find_false_negatives()
        
        self.calculate_sensitivity()

        self.calculate_precision()
        
        return true_positives_cluster_file_name


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
    category_file.write("#cluster_type	gene1	gene2	max_split_cnt	max_span_cnt	sample_type	disease	tools	inferred_fusion_type	gene1_on_bnd	gene1_close_to_bnd	gene2_on_bnd	gene2_close_to_bnd	dna_supp	samples	chr1	breakpoint_1	chr2	breakpoint_2	gene1_candidates	gene2_candidates	gene1_strands	gene2_strands"  + "\n")
    for category_fusion in filtered_list:
        category_file.write(category_fusion.line + "\n")
    category_file.close()

    return output_file





#**********************
# PERFORM VALIDATION
#**********************

#cluster_files output directory
cluster_outdir = args.output_dir + "/" + "cluster_stats_files" + "/" 

# create ValidatedFusionStats object for overall .cluster file 
total_fusion_stats = ValidatedFusionStats(cluster_outdir, args.validated_fusions_file, args.fusion_cluster_file, "total")

#TEST
print(total_fusion_stats.samples)
print(total_fusion_stats.num_samples)
#TEST

#validation summary output file
#stats_output_file_name = args.output_dir + "/" + "validation_output_stats." + args.sample_name + ".txt"
stats_output_file_name = args.output_dir + "/" + "validation_output_stats.txt"
# erase old data from this file by writing a blank line
stats_output_file = open(stats_output_file_name, 'w+')
stats_output_file.write("")
stats_output_file.close()


# Divide cluster file into 4 separate files create new ValidatedFusionStats object for each tool
fusion_stats_objects = []

all_cluster_files = []
# generate .cluster files for each tool
#tools = ["defuse", "integrate", "fusionmap", "ericscript", "STAR-Fusion"]
tools = ["defuse", "integrate", "fusionmap", "ericscript", "star_fusion"]
#tools = ["defuse", "integrate", "fusionmap", "ericscript"]
for tool in tools:
    #generate category file
    category_file = generate_filtered_category_file(total_fusion_stats, cluster_outdir + tool + ".cluster", tool=tool)
    all_cluster_files.append(cluster_outdir + tool + ".cluster")
    #create fusion stats object
    fusion_stats_object = ValidatedFusionStats(cluster_outdir, args.validated_fusions_file, category_file, tool)
    fusion_stats_objects.append(fusion_stats_object)

# TWO OR MORE
# generate .cluster file for 2 or more tools
two_or_more_category_file = generate_filtered_category_file(total_fusion_stats, cluster_outdir + "merged.two_or_more"+ ".cluster", num=2)
two_or_more_fusion_stats = ValidatedFusionStats(cluster_outdir, args.validated_fusions_file, two_or_more_category_file, "two_or_more")
fusion_stats_objects.append(two_or_more_fusion_stats)

# LOOK AT TWO OR MORE FOR EACH TOOL SEPARATELY
for tool in tools:
    category_file = generate_filtered_category_file(two_or_more_fusion_stats, cluster_outdir + tool + "two_or_more.cluster", tool=tool)
    fusion_stats_object = ValidatedFusionStats(cluster_outdir, args.validated_fusions_file, category_file, tool + "-two_or_more")
    fusion_stats_objects.append(fusion_stats_object)


# THREE OR MORE
# generate .cluster file for 3 or more tools
three_or_more_category_file = generate_filtered_category_file(total_fusion_stats, cluster_outdir + "merged.three_or_more"+ ".cluster", num=3)
three_or_more_fusion_stats = ValidatedFusionStats(cluster_outdir, args.validated_fusions_file, three_or_more_category_file, "three_or_more")
fusion_stats_objects.append(three_or_more_fusion_stats)

#STAR-Fusion + two more 
#STAR_Fusion_plus_two_category_file = generate_filtered_category_file(three_or_more_fusion_stats, cluster_outdir + "STAR-Fusion_plus_two" + ".cluster", tool="STAR-Fusion")
STAR_Fusion_plus_two_category_file = generate_filtered_category_file(three_or_more_fusion_stats, cluster_outdir + "STAR-Fusion_plus_two" + ".cluster", tool="star_fusion")
STAR_Fusion_plus_two_fusion_stats = ValidatedFusionStats(cluster_outdir, args.validated_fusions_file, STAR_Fusion_plus_two_category_file, "STAR-Fusion_plus_two")
fusion_stats_objects.append(STAR_Fusion_plus_two_fusion_stats)

# FOUR OR MORE
# generate .cluster file for 4 or more tools
four_or_more_category_file = generate_filtered_category_file(total_fusion_stats, cluster_outdir + "merged.four_or_more"+ ".cluster", num=4)
four_or_more_fusion_stats = ValidatedFusionStats(cluster_outdir, args.validated_fusions_file, four_or_more_category_file, "four_or_more")
fusion_stats_objects.append(four_or_more_fusion_stats)

# FIVE OR MORE
# generate .cluster file for 4 or more tools
five_or_more_category_file = generate_filtered_category_file(total_fusion_stats, cluster_outdir + "merged.five_or_more"+ ".cluster", num=5)
five_or_more_fusion_stats = ValidatedFusionStats(cluster_outdir, args.validated_fusions_file, five_or_more_category_file, "five_or_more")
fusion_stats_objects.append(five_or_more_fusion_stats)

# TOTAL
fusion_stats_objects.append(total_fusion_stats)

# generate stats
stats_output_file = open(stats_output_file_name, 'a')
true_positives_cluster_files = []
stats_output_file.write('Tool	num_detected_fusions	num_TP	num_unvalidated	num_FN	sensitivity	precision\n')
for fusion_stats in fusion_stats_objects:
    #print "generating cluster files for " + fusion_stats_object.name
    # generate TP, FP and FN files
    true_positives_cluster_file = fusion_stats.compare_validated_and_output_fusions()
    true_positives_cluster_files.append(true_positives_cluster_file)
    #print '{}	{}	{}'.format(fusion_stats.name, fusion_stats.num_fusions, fusion_stats.num_true_positives) 
    stats_output_file.write('{}	{}	{}	{}	{}	{}	{}\n'.format(fusion_stats.name, fusion_stats.num_fusions, fusion_stats.num_true_positives, fusion_stats.num_unvalidated_fusions, fusion_stats.num_false_negatives, fusion_stats.sensitivity, fusion_stats.precision) )
stats_output_file.close()

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

# GENERATE TRUE POSITIVE CATEGORY COUNTS FOR EACH CALLER
#fusionmap.true_positives.cluster
category_count_file = open(os.path.join(args.output_dir, "true_positive.category_count_file.txt"), 'w+')
categories = ['GeneFusion', 'SameGene', 'NoDriverGene', 'ReadThrough', 'TruncatedNoncoding', 'TruncatedCoding']
category_count_file.write("\t".join(["tool"] + categories) + "\n" )

for tool in tools:
    fusion_stats = CategoryFusionStats(os.path.join(cluster_outdir, tool + ".true_positives.cluster"))
    category_dict = fusion_stats.generate_category_counts()
    print(tool, "Number of fusions:", fusion_stats.num_fusions, category_dict)
    category_count_file.write("\t".join([tool] + [str(category_dict[category]) for category in categories] ) + "\n" )
category_count_file.close()



# ***********
# SCRAP BELOW
# ***********

# SUMMARY FILE
#def generate_summary_file(merged_TP_list ):
#    """
#    generate summary file which contains only the folowing columns: gene1   gene2   max_split_cnt   max_span_cnt    tools   samples chr1    breakpoint_1    chr2    breakpoint_2
#    """
#    merged_TP_summary_file_name = "output/"+ str(merged_TP_list[0].samples[0]) + ".total.true_positives.summary"
#    #print(merged_TP_summary_file_name)
#    merged_TP_summary_file = open(merged_TP_summary_file_name, 'w+')
#    merged_TP_summary_file.write("gene1	gene2	max_split_cnt	max_span_cnt	tools	samples	chr1	breakpoint_1	chr2	breakpoint_2\n")	
#    for fusion in merged_TP_list:
#        # output line of category_fusion using class attributes, since original category_fusion.line string is no longer accurate
#        merged_TP_summary_file.write("\t".join(map(str, [fusion.gene1, fusion.gene2, fusion.max_split_cnt, fusion.max_span_cnt, ",".join(list(set(sorted(fusion.tools)))), ",".join(list(set(fusion.samples))), fusion.chr1, ",".join(list(set([str(bp) for bp in fusion.breakpoint_1]))), fusion.chr2, ",".join(list(set([str(bp) for bp in fusion.breakpoint_2])))])) + "\n")
#
#generate_summary_file(merged_TP_list)



# Generate true_positive .cluster file merged based on gene names
#merged_TP_cluster_file_name = "cluster_outdir/merged.true_positives.cluster"
#merged_TP_list = generate_merged_category_list_by_gene_names(true_positives_cluster_files)
#generate_merged_TP_cluster_file()



# generate ValidatedFusionStats objects for both merged_TP and merged_total
#merged_TP_fusion_stats = ValidatedFusionStats(args.validated_fusions_file, merged_TP_cluster_file_name, "merged_TP_cluster") 
#merged_total_fusion_stats = ValidatedFusionStats(args.validated_fusions_file, merged_total_cluster_file_name, "merged_total_cluster") 
#print >> sys.stderr, "merged_TP_fusions" + str([fusion.line for fusion in merged_TP_fusion_stats.category_list])

# TWO OR MORE

# generate TP and total .cluster files for 2 or more tools
#two_or_more_TP_category_file = generate_filtered_category_file(merged_TP_fusion_stats, "cluster_outdir/" + "merged.two_or_more_TP" + ".cluster", num=2)
#two_or_more_total_category_file = generate_filtered_category_file(merged_total_fusion_stats, "cluster_outdir/" + "merged.two_or_more_total" + ".cluster", num=2)

# generate ValidatedFusionStats objects for both merged.two_or_more_TP and merged.two_or_more_total 
#two_or_more_TP_fusion_stats = ValidatedFusionStats(args.validated_fusions_file, two_or_more_TP_category_file, "two_or_more")
#two_or_more_total_fusion_stats = ValidatedFusionStats(args.validated_fusions_file, two_or_more_total_category_file, "two_or_more")


# print line of output file

#stats_output_file = open(stats_output_file_name, 'a')
#generate_two_or_more_TP_fusion_stats(two_or_more_TP_fusion_stats, two_or_more_total_fusion_stats)
#stats_output_file.close()
## TOTAL
#stats_output_file = open(stats_output_file_name, 'a')
#generate_merged_TP_fusion_stats(merged_TP_fusion_stats, merged_total_fusion_stats)
#stats_output_file.close()






# generate file containing true positive fusions detected by 2 or more tools in validated_fusion_file format
#two_or_more_valfile(fusion_stats_objects)


#def two_or_more_valfile(fusion_stats_objects):
#    """
#    generate file containing true positive fusions detected by 2 or more tools in validated_fusion_file format
#    Uses the unique string present in the validation file. Does not use gene names
#    """
#    # combine true_positive_validated_fusions lists from each tool into one giant list
#    sublists = [fusion_stats.true_positive_validated_fusions for fusion_stats in fusion_stats_objects]
#    # flatten sublists
#    master_list =  [item.line for sublist in sublists for item in sublist]
#    
#    #master_list = [val_fusion.line for val_fusion in fusion_object_list ]
#    
#    # count the number of times each line is present in master_list
#    count_dict = {}
#    for i in master_list:
#        count_dict[i] = count_dict.get(i, 0) + 1
#    
#    # return sublist of lines that are present 2 or more times
#    two_or_more = [line for line in count_dict.keys() if count_dict[line] >1]
#    print("two_or_more	NaN	{}".format(len(two_or_more)))
#    print("all 	NaN	{}".format(len(count_dict.keys())))

#def generate_two_or_more_TP_fusion_stats(TP, total):
#    """
#    Print stats line for two_or_more_TP
#    """
#    # run validation for both objects
#    temp = TP.compare_validated_and_output_fusions()
#    temp = total.compare_validated_and_output_fusions()
#    # print hybrid line for "two_or_more"
#    num_unvalidated_fusions = total.num_fusions - TP.num_true_positives
#    sensitivity = float(TP.num_true_positives)/float(len(TP.validated_fusions))
#    try:
#        precision = float(TP.num_true_positives)/float(TP.num_true_positives + num_unvalidated_fusions)
#    except ZeroDivisionError:
#        precision = 0
#    stats_output_file.write('{} {}      {}      {}      {}      {}      {}\n'.format(TP.name, total.num_fusions, TP.num_true_positives, num_unvalidated_fusions, TP.num_false_negatives, sensitivity, precision) )
#
#def generate_merged_TP_fusion_stats(TP, total):
#    # run validation for both objects
#    temp = TP.compare_validated_and_output_fusions()
#    temp = total.compare_validated_and_output_fusions()
#
#    # print line of output file
#    num_unvalidated_fusions = total.num_fusions - TP.num_true_positives
#    sensitivity = float(TP.num_true_positives)/float(len(TP.validated_fusions))
#    precision = float(TP.num_true_positives)/float(TP.num_true_positives + num_unvalidated_fusions)
#    stats_output_file.write('{}       {}      {}      {}      {}      {}      {}\n'.format("total", total.num_fusions, TP.num_true_positives, num_unvalidated_fusions, TP.num_false_negatives, sensitivity, precision) ) 

#def map_fusion_genes_to_bed_file(fusion):
#    """
#    .cluster file fields:
#    cluster_type    gene1   gene2  max_split_cnt max_span_cnt sample_type disease tools inferred_fusion_type gene1_on_bnd gene1_close_to_bnd  gene2_on_bnd gene2_close_to_bnd dna_supp  samples chr1 breakpoint_1 chr2 breakpoint_2 
#    """
#    left = []
#    right = []
#    for genebed in genebed_objects:
#        if ( ( genebed.start <= fusion.breakpoint_1 <= genebed.end)) and fusion.chr1 == genebed.chr:
#            fusion_left_names = [tup[0] for tup in left]
#            if genebed.gene_name not in fusion_left_names:
#                left.append((genebed.gene_name, genebed.type, genebed.strand, [genebed.start, genebed.end]))
#                left.append((genebed.gene_name, genebed.type, genebed.strand, [genebed.start, genebed.end]))
#        if ( ( genebed.start <= fusion.breakpoint_2 <= genebed.end)) and fusion.chr2 == genebed.chr:
#            fusion_left_names = [tup[0] for tup in left]
#            if genebed.gene_name not in fusion_left_names:
#                left.append((genebed.gene_name, genebed.type, genebed.strand, [genebed.start, genebed.end]))
#

#def generate_merged_TP_cluster_file():
#    # Generate true_positive .cluster file merged based on gene names
#    merged_TP_cluster_file = open(merged_TP_cluster_file_name, 'w+')
#    merged_TP_cluster_file.write("#cluster_type    gene1   gene2  max_split_cnt max_span_cnt sample_type disease tools inferred_fusion_type gene1_on_bnd gene1_close_to_bnd  gene2_on_bnd gene2_close_to_bnd dna_supp  samples chr1 breakpoint_1 chr2 breakpoint_2"  + "\n")
#    for fusion in merged_TP_list:
#        # output line of category_fusion using class attributes, since original category_fusion.line string is no longer accurate
#         merged_TP_cluster_file.write("\t".join(map(str, [fusion.cluster_type, fusion.gene1, fusion.gene2, fusion.max_split_cnt, fusion.max_span_cnt, fusion.sample_type, ",".join(list(set(fusion.disease))), ",".join(list(set(sorted(fusion.tools)))), fusion.inferred_fusion_type, fusion.gene1_on_bnd, fusion.gene1_close_to_bnd, fusion.gene2_on_bnd, fusion.gene2_close_to_bnd, fusion.dna_supp, ",".join(list(set(fusion.samples))), fusion.chr1, ",".join(list(set([str(bp) for bp in fusion.breakpoint_1]))), fusion.chr2, ",".join(list(set([str(bp) for bp in fusion.breakpoint_2]))), fusion.gene1_candidates, fusion.gene2_candidates, fusion.gene1_strands, fusion.gene2_strands])) + "\n")
#    merged_TP_cluster_file.close()
