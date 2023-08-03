#!/usr/bin/env python


debug=0
import sys
import pysam
import bisect
import time
import sequtils
import copy
import re
#import pybedtools.bedtool


# fusions are clusted by gene pairs
# Category file format:
#Gene_Cluster    HES1    IFNL1   43      41      Tumor   VALIDATION      integrate,defuse,ericscript     GeneFusion      True    True    True    True    -1      smc_rna_sim45   chr3    193854837       chr19   39787445

class CategoryFusions():
    """
    Represents a line of the output file of the fusion pipeline. Each line of the file represents a possible fusion event.
    The output file is named merged.cff.reann.dnasupp.bwafilter.30.cluster
    """
    def __init__(self, category_line):
        self.__load_category(category_line)
    def __load_category(self, category_line):
        tmp = category_line.split()
        #HEADER
        #cluster_type gene1 gene2 max_split_cnt max_span_cnt sample_type disease tools inferred_fusion_type gene1_on_bnd gene1_close_to_bnd gene2_on_bnd gene2_close_to_bnd dna_supp samples chr1 breakpoint_1 chr2 breakpoint_2 captured_reads_tumor_mean captured_reads_normal_mean fusion_IDs
        self.cluster_type = tmp[0]
        self.gene1 = tmp[1]
        self.gene2 = tmp[2]
        self.max_split_cnt = int(tmp[3])
        self.max_span_cnt = int(tmp[4])
        self.sample_type = tmp[5] 
        self.disease = tmp[6].split(",")
        #self.tools = "ericscript")
        self.tools = tmp[7].split(",")
        self.inferred_fusion_type = tmp[8]
        self.gene1_on_bnd = tmp[9]
        self.gene1_close_to_bnd = tmp[10]
        self.gene2_on_bnd = tmp[11]
        self.gene2_close_to_bnd = tmp[12]   
        self.cancer_db_hits = tmp[13].split(',')
        # dna support values: 
        #0: has data, no dna pairs; >0:  number of dna pair clusters; -1: no data;
        # -2: gene not in annotation; -3: chr not in bam; -4: confilicting window start and end 
        self.samples = tmp[14].split(",")
        self.line = category_line.strip()
        if len(tmp) > 15:
            # add 4 new columns for purposes of validation of chromosomal positions with reference to validated fusions
            self.chr1 = tmp[15]
            #setting breakpoint lists 
            self.breakpoint_1 = [int(num) for num in tmp[16].split('|')]
            self.chr2 = tmp[17]
            self.breakpoint_2 = [int(num) for num in tmp[18].split('|')]
            self.tools_breakpoints_dict_1 = dict( zip([tool for tool in self.tools], [breakpoint for breakpoint in self.breakpoint_1 ]) )
            self.tools_breakpoints_dict_2 = dict( zip([tool for tool in self.tools], [breakpoint for breakpoint in self.breakpoint_2 ]) )
        # attributes to accrue list of gene names in bed feature file that intersect with breakpoints 
        self.left = []
        self.right = []
        self.exon1 = tmp[19]
        self.exon2 = tmp[20]
        self.fusion_IDs = tmp[21].split(',')
       # try:
       #     self.gene1_candidates = tmp[19] 
       # except IndexError:
       #     self.gene1_candidates = ""
       # try:
       #     self.gene2_candidates = tmp[20]
       # except IndexError:
       #     self.gene2_candidates = ""
       # try:
       #     # GENE STRANDS REPRESENT CANDIDATE STRANDS QUERIED FROM GENE ANNOTATION FILE
       #     self.gene1_strands = tmp[21]
       # except IndexError:
       #     self.gene1_strands = ""
       # try:
       #     self.gene2_strands = tmp[22]
       # except IndexError:
       #     self.gene2_strands = "" 

    def out(self):
        #print self.line
        print "\t".join(map(str, [self.cluster_type, self.gene1, self.gene2, self.max_split_cnt, self.max_span_cnt, self.sample_type, ",".join(self.disease), ",".join(self.tools), self.inferred_fusion_type, self.gene1_on_bnd, self.gene1_close_to_bnd, self.gene2_on_bnd, self.gene2_close_to_bnd, ",".join(self.cancer_db_hits), ",".join(self.samples), self.chr1, "|".join(map(str, self.breakpoint_1)), self.chr2, "|".join(map(str, self.breakpoint_2)), self.exon1, self.exon2, ",".join(self.fusion_IDs)]))

    def out_subset(self):
        # gene1   gene2   chr1    breakpoint_1    chr2    breakpoint_2    max_split_cnt   max_span_cnt    sample_type disease tools   inferred_fusion_type    samples fusion_IDs
        print "\t".join(map(str, [self.gene1, self.gene2,  self.chr1, "|".join(map(str, self.breakpoint_1)), self.chr2, "|".join(map(str, self.breakpoint_2)), self.max_split_cnt, self.max_span_cnt, self.sample_type, ",".join(self.disease), ",".join(self.tools), self.inferred_fusion_type, ",".join(self.samples), ",".join(self.cancer_db_hits),",".join(self.fusion_IDs)]))

class CategoryFusionSubset():
    """
    A subset of the .cluster file, abridged for easy reading
    """
    def __init__(self, category_line):
        self.__load_category(category_line)
    def __load_category(self, category_line):
        tmp = category_line.split()
        self.line = category_line.strip()
        #HEADER
        # gene1   gene2   chr1    breakpoint_1    chr2    breakpoint_2    max_split_cnt   max_span_cnt    sample_type disease tools   inferred_fusion_type    samples cancer_db_hits captured_reads_normal_mean fusion_IDs
        self.gene1 = tmp[0]
        self.gene2 = tmp[1]
        self.chr1 = tmp[2]
        self.breakpoint_1 = [int(num) for num in tmp[3].split('|')]
        self.chr2 = tmp[4]
        self.breakpoint_2 = [int(num) for num in tmp[5].split('|')]
        self.max_split_cnt = int(tmp[6])
        self.max_span_cnt = int(tmp[7])
        self.sample_type = tmp[8] 
        self.disease = tmp[9].split(",")
        self.tools = tmp[10].split(",")
        self.inferred_fusion_type = tmp[11]
        self.samples = tmp[12].split(",")
        self.cancer_db_hits = tmp[13].split(',') 
        self.fusion_IDs = tmp[14].split(',')

    def out_subset(self):
        # gene1   gene2   chr1    breakpoint_1    chr2    breakpoint_2    max_split_cnt   max_span_cnt    sample_type disease tools   inferred_fusion_type    samples cancer_db_hits captured_reads_normal_mean fusion_IDs
        print "\t".join(map(str, [self.gene1, self.gene2,  self.chr1, "|".join(map(str, self.breakpoint_1)), self.chr2, "|".join(map(str, self.breakpoint_2)), self.max_split_cnt, self.max_span_cnt, self.sample_type, ",".join(self.disease), ",".join(self.tools), self.inferred_fusion_type, ",".join(self.samples), ",".join(self.cancer_db_hits), ",".join(self.fusion_IDs)]))


#HEAD OUTPUT FUNCTIONS
def output_cluster_header():
        print "\t".join(["#cluster_type", "gene1", "gene2", "max_split_cnt", "max_span_cnt", "sample_type", "disease", "tools", "inferred_fusion_type", "gene1_on_bnd", "gene1_close_to_bnd", "gene2_on_bnd", "gene2_close_to_bnd", "cancer_db_hits", "samples", "chr1", "breakpoint_1", "chr2", "breakpoint_2", "exon1", "exon2", "fusion_IDs"])
def output_cluster_header_subset():
        print "\t".join(["#gene1", "gene2", "chr1", "breakpoint_1", "chr2", "breakpoint_2", "max_split_cnt", "max_span_cnt", "sample_type", "disease", "tools", "inferred_fusion_type", "samples", "cancer_db_hits", "fusion_IDs"])
        
            
class CategoryFusionStats():
    """
    Represents an entire category file, the output file of the fusion pipeline.
    Has an instance variable, self.category_list, which is a list of CategoryFusion objects. Each CagegoryFusion
    object represents a line of the category file.
    Contains various methods for filtering the category_file, each which returns a subset of the fusion
    events contained in the category file.
    """
    def __init__(self, category_file):
        self.category_list = []
        self.num_fusions = 0
        #Check if Subset
        file = open(category_file, "r")
        line = file.readline().split("\t")
        file.close()
        if line[0] == "#gene1": self.__load_category_file_subset(category_file) 
        else: self.__load_category_file(category_file)

    def __load_category_file(self, category_file):
        for line in open(category_file, "r"):
            if line.startswith("#"):
                continue
            category = CategoryFusions(line)
            self.category_list.append(category)
        self.num_fusions = len(self.category_list)
        # count number of samples present in category file
        samples = []
        for fusion in self.category_list:
            samples += fusion.samples
        self.samples = set(samples)
        self.num_samples = len(self.samples)

    def __load_category_file_subset(self, category_file):
        for line in open(category_file, "r"):
            if line.startswith("#"):
                continue
            category = CategoryFusionSubset(line)
            self.category_list.append(category)
        self.num_fusions = len(self.category_list)
        # count number of samples present in category file
        samples = []
        for fusion in self.category_list:
            samples += fusion.samples
        self.samples = set(samples)
        self.num_samples = len(self.samples)

    def filter_clustered_fusions(self, fusion_list, split=0, span=0, tool_num=0, on_bnd=False, close_to_bnd=False, fusion_type=[]):
        filtered_list = []
        for fu in fusion_list:
            keep = False
            if fu.max_split_cnt >= split and fu.max_span_cnt >= span and len(fu.tools) >= tool_num:
                if on_bnd:
                    if fu.gene1_on_bnd and fu.gene2_on_bnd:
                        keep = True
                elif close_to_bnd:
                    if fu.gene1_close_to_bnd and fu.gene2_close_to_bnd:
                        keep = True
                else:
                        keep = True
            if keep:
                filtered_list.append(fu)
        return filtered_list
                
    def filter_samples(self, fusion_list, sample_list):
        return [x for x in fusion_list if set(sample_list).intersection(set(x.samples))]

    def filter_recurrent(self, fusion_list, threshold):
        return filter(lambda x:len(x.samples)>=threshold, fusion_list)

    def filter_disease(self, fusion_list, disease):
        return filter(lambda x:disease in x.disease, fusion_list)

    def filter_sample_number(self, fusion_list, threshold, sample_prefix):
        return filter(lambda x:len(filter(lambda y:y.startswith(sample_prefix), x.samples)) >= threshold, fusion_list)

    def filter_split_cnt (self, fusion_list, split_cnt):    
        return filter(lambda x:x.max_split_cnt>=split_cnt, fusion_list)

    def filter_span_cnt (self, fusion_list, span_cnt):  
        return filter(lambda x:x.max_span_cnt>=span_cnt, fusion_list)

    def filter_sample_type (self, fusion_list, sample_type):    
        return filter(lambda x:x.sample_type==sample_type, fusion_list)

    def filter_tools_name (self, fusion_list, tool_name):   
        return filter(lambda x:tool_name in x.tools, fusion_list)
    
    def filter_tools_name_multi (self, fusion_list, tool_names):  
         
        return filter(lambda x:sorted(tool_names) == sorted(x.tools), fusion_list)

    def filter_tools_num (self, fusion_list, tools_num):    
        return filter(lambda x:len(x.tools)>=tools_num, fusion_list)

    def filter_inferred_type (self, fusion_list, inferred_type):    
        return filter(lambda x:x.inferred_fusion_type==inferred_type, fusion_list)

    def filter_on_bnd (self, fusion_list):  
        return filter(lambda x:x.gene1_on_bnd and x.gene2_on_bnd, fusion_list)

    def filter_close_to_bnd (self, fusion_list):    
        return filter(lambda x:(x.gene1_close_to_bnd == "True" and x.gene2_close_to_bnd == "True"), fusion_list)

    def filter_dna_supp (self, fusion_list):    
        return filter(lambda x:int(x.dna_supp)>0, fusion_list)

    def generate_category_counts(self):
        """
        Counts fusions for each category (i.e. self.inferred_fusion_type )
        return: dictionary of { CategoryName: num_fusions }
        """
        category_dict = {}
        categories = ['SameGene', 'NoDriverGene', 'ReadThrough', 'GeneFusion', 'TruncatedNoncoding', 'TruncatedCoding']
        for category in categories:
            category_dict[category] = 0 
        for fusion in self.category_list:
            if category_dict[fusion.inferred_fusion_type] == 0:
                category_dict[fusion.inferred_fusion_type] = 1
            else: 
                category_dict[fusion.inferred_fusion_type] += 1          
        return category_dict

class CffFusionStats():
    __fusion_dict = {}
    __fusion_samples_dict = {}
    def __init__(self, cff_file):
        pass
    def __load_fusions(self, cff_file):
        for line in open(cff_file, "r"):
            fusion = CffFusion(line)
            key = fusion.t_gene1 + "_" + fusion.t_gene2
            rkey = fusion.t_gene2 + "_" + fusion.t_gene1
            if rkey not in self.__fusion_dict:
                self.__fusion_dict.setdefault(key, []).append(fusion)
                self.__fusion_samples_dict.setdefault(key, []).append(fusion.sample_name)
            else:
                self.__fusion_dict.setdefault(rkey,[]).append(fusion)
                self.__fusion_samples_dict.setdefault(rkey, []).append(fusion.sample_name)
    
    # compare two fusions based on the up/downstream gene sets, if overlap on both sets consider them as same fusion
    def is_same_gene_pair_fusion(self, fusion1, fusion2):
        id = 1 
        # compare gene sets on fw strand 
        up_g1, down_g1 = set(fusion1.get_reannotated_genes(id))
        up_g2, down_g2 = set(fusion2.get_reannotated_genes(id))

        if up_g1 & up_g2 and down_g1 & down_g2:
            return True

        id = 2 
        # compare gene sets on bw strand 
        up_g1, down_g1 = set(fusion1.get_reannotated_genes(id))
        up_g2, down_g2 = set(fusion2.get_reannotated_genes(id))

        if up_g1 & up_g2 and down_g1 & down_g2:
            return True

        return False

    # compare gene fusions based on breakpoints
    def generate_common_fusion_stats(self, cff_file):
        fusion_bp_dict = {}
        fusion_pos1_idx_dict = {}
        # build index for fusions based on fusion.pos1
        for line in  open(cff_file, "r"):
            if line.startswith("#"):
                continue
            fusion = CffFusion(line)
            fusion_bp_dict.setdefault(fusion.chr1, {}).setdefault(fusion.strand1, []).append(fusion)
        for chr in fusion_bp_dict:
            for strand in fusion_bp_dict[chr]:
                fusion_bp_dict[chr][strand].sort(key = lambda x:x.pos1) 
                fusion_pos1_idx_dict.setdefault(chr, {}).setdefault(strand, [f.pos1 for f in fusion_bp_dict[chr][strand]]) 
        # find common fusions

    # output fusions in a fusion list as clustered fusions
    def output_clustered_fusions(self, fusion_list, cluster_type):
            use_reference_fusion=1
            #use first fusion in fusion as reference fusion to account for flipped fusions for BP_Cluster
            if use_reference_fusion:
                ref_fus = fusion_list[0]
                #ref_fus_bp1 = (ref_fus.chr1, ref_fus.pos1, ref_fus.strand1)
                #flip all fusions not ordered the same way as reference:
                for fusion in fusion_list:
                    #TEST: remove reannotated gene name
                    #fusion.reann_gene1 = fusion.t_gene1
                    #fusion.reann_gene2 = fusion.t_gene2

                    #fus_bp1 = (fusion.chr1, fusion.pos1, fusion.strand1)
                    # check both breakpoint distance AND gene names
                    #if not cmp_fusion_breakpoints(ref_fus_bp1, fus_bp1, 100000) and not (ref_fus.reann_gene1 == fusion.reann_gene1):
                    if not (ref_fus.t_gene1 == fusion.t_gene1):
                        #flip fusion:
                        fusion.chr1, fusion.chr2 = fusion.chr2, fusion.chr1
                        fusion.pos1, fusion.pos2 = fusion.pos2, fusion.pos1
                        fusion.strand1, fusion.strand2 = fusion.strand2, fusion.strand1
                        fusion.t_gene1, fusion.t_gene2 = fusion.t_gene2, fusion.t_gene1
                        fusion.t_area1, fusion.t_area2 = fusion.t_area2, fusion.t_area1
                        #fusion.reann_gene1, fusion.reann_gene2 = fusion.reann_gene2, fusion.reann_gene1
                        fusion.reann_type1,fusion.reann_type2 = fusion.reann_type2,fusion.reann_type1

            #need to remove NA values from below lists, but only if there are valid gene names in them
            #gene1_list_orig = [f.reann_gene1 for f in fusion_list]
            gene1_list_orig = [f.t_gene1 for f in fusion_list]
            gene1_list = [gene for gene in gene1_list_orig if gene !="NA"]
            if len(gene1_list) < 1:
                gene1_list = gene1_list_orig
            #gene2_list_orig = [f.reann_gene2 for f in fusion_list]
            gene2_list_orig = [f.t_gene2 for f in fusion_list]
            gene2_list = [gene for gene in gene2_list_orig if gene !="NA"]
            if len(gene2_list) < 1:
                gene2_list = gene2_list_orig

            max_split_cnt = max([f.split_cnt for f in fusion_list])
            max_span_cnt = max([f.span_cnt for f in fusion_list])
            sample_list = [f.sample_name for f in fusion_list]  

            # get fusion.captured_reads averages for T and N
            captured_reads_tumor = [f.captured_reads for f in fusion_list if f.sample_type == "Tumor"] 
            try:
                captured_reads_tumor_mean = sum(captured_reads_tumor)/float(len(captured_reads_tumor))
            except ZeroDivisionError:
                captured_reads_tumor_mean = -1
            captured_reads_normal = [f.captured_reads for f in fusion_list if f.sample_type == "Normal"] 
            try:
                captured_reads_normal_mean = sum(captured_reads_normal)/float(len(captured_reads_normal))
            except ZeroDivisionError:
                captured_reads_normal_mean = -1

            # list of Fusion IDs
            fusion_IDs = [f.fusion_id for f in fusion_list] 
            disease_list = [f.disease for f in fusion_list] 
            tool_list = [f.tool for f in fusion_list]   
            sample_type_list = [f.sample_type for f in fusion_list] 
            gene1_on_bndry = "True" in [f.gene1_on_bndry for f in fusion_list] 
            gene1_close_to_bndry = "True" in [f.gene1_close_to_bndry for f in fusion_list] 
            gene2_on_bndry = "True" in [f.gene2_on_bndry for f in fusion_list] 
            gene2_close_to_bndry = "True" in [f.gene2_close_to_bndry for f in fusion_list] 
            
            dna_supp_cluster_num = max([int(f.dnasupp) for f in fusion_list])
            
            category_list = [f.category for f in fusion_list]   
            # PRIORITIZE CATEGORIES TO REMOVE MULTIPLE CATEGORIES
            if "ReadThrough" in category_list:
                category_list = ["ReadThrough"]
            elif "CodingFusion" in category_list:
                category_list = ["CodingFusion"]
            elif "TruncatedCoding" in category_list:
                category_list = ["TruncatedCoding"]
            elif "TruncatedNoncoding" in category_list:
                category_list = ["TruncatedNoncoding"]
            elif "NoDriverGene" in category_list:
                category_list = ["NoDriverGene"]
            
            # added 4 new fields to final .cluster file for purposes of validation 
            chr1_list = [str(f.chr1) for f in fusion_list]
            breakpoint_1_list = [str(f.pos1) for f in fusion_list]
            chr2_list = [str(f.chr2) for f in fusion_list]
            breakpoint_2_list = [str(f.pos2) for f in fusion_list]
            # print statement modified to include the 4 above new fields
            print "\t".join(map(str, [cluster_type, ",".join(list(set(gene1_list))), ",".join(list(set(gene2_list))), max_split_cnt, max_span_cnt, ",".join(list(set(sample_type_list))), ",".join(list(set(disease_list))), ",".join(list(set(tool_list))), ",".join(list(set(category_list))), gene1_on_bndry, gene1_close_to_bndry, gene2_on_bndry, gene2_close_to_bndry, dna_supp_cluster_num, ",".join(list(set(sample_list))), ",".join(list(set(chr1_list))), ",".join(list(set(breakpoint_1_list))), ",".join(list(set(chr2_list))), ",".join(list(set(breakpoint_2_list))), captured_reads_tumor_mean, captured_reads_normal_mean,",".join(list(set(fusion_IDs)))]))


    # cluster fusions of "NoDriverGene" and "Truncated" type on their breakpoints
    def generate_common_fusion_stats_by_breakpoints(self, fusion_list):
        diff = 100000
        # save clustered fusion id, skip a fusion when it is already clustered
        clustered_id = {}
        print >> sys.stderr, "bp fusion list:", len(fusion_list)
        for i in range(len(fusion_list)):
            if i in clustered_id:
                continue

            fusion1 = fusion_list[i]
            #ORDERED BREAKPOINTS ALLOWS FUSIONS WITH SAME GENE PAIRS TO BE ADDED TO gene_cluster_list,
            # BUT FUNCTION "output_clustered_fusions" does not account for the fact that gene orders might
            # be flipped when outputting a cluster
            small_bp1, big_bp1 = fusion1.get_ordered_breakpoints()
            clustered_id.setdefault(i, i)
            fusion_cluster_list = []
            fusion_cluster_list.append(fusion1)
            for j in range(len(fusion_list)):
                if j in clustered_id:
                    continue
                if i == j:
                    continue

                fusion2 = fusion_list[j]
                small_bp2, big_bp2 = fusion2.get_ordered_breakpoints()
                if cmp_fusion_breakpoints(small_bp1, small_bp2, diff) and cmp_fusion_breakpoints(big_bp1, big_bp2, diff):
                    clustered_id.setdefault(j, j)
                    fusion_cluster_list.append(fusion2)

            self.output_clustered_fusions(fusion_cluster_list, "BP_Cluster")
                

    # compare gene fusions based on re-annotated gene names, hard to work for truncated fusions, used generate_common_fusion_stats_by_breakpoints for NoDriverGenes and Trucated type fusions
    def generate_common_fusion_stats_by_genes(self, cff_file):
        fusion_dict = {}
        fusion_list_for_bp_cmp = []
        common_key_dict = {}
        # cluster fusions by gene pairs, save in fusion_dict
        for line in open(cff_file, "r"):
            if line.startswith("#"):
                continue
            fusion = CffFusion(line)
            # send fusion to breakpoint cluster later
            #if fusion.reann_gene1 == "NA" or fusion.reann_gene2 == "NA":
            if fusion.t_gene1 == "NA" or fusion.t_gene2 == "NA":
                fusion_list_for_bp_cmp.append(fusion)
            else:
                #key = (fusion.reann_gene1 + "|" + fusion.chr1, fusion.reann_gene2+ "|" + fusion.chr2)  
                key = ",".join(sorted([fusion.t_gene1 + "|" + fusion.chr1, fusion.t_gene2+ "|" + fusion.chr2])) 
                fusion_dict.setdefault(key, []).append(fusion)
        # output clustered fusions
        for key in fusion_dict:
            fusion_list = fusion_dict[key]

            self.output_clustered_fusions(fusion_list, "Gene_Cluster")
        self.generate_common_fusion_stats_by_breakpoints(fusion_list_for_bp_cmp)        
        
    def get_gene_order_stats(self):
        for key in self.__fusion_dict:
            n_sg = 0 # same gene
            n_rt = 0 # read through
            n_gf = 0 # gene fusion
            n_tc = 0 # truncated coding
            n_tn = 0 # truncated noncoding
            n_ns = 0 # nonsense
            gene_order = []
            sample_type = []
            for sample in set(self.__fusion_samples_dict[key]):
                if sample.endswith("N"):
                    sample_type.append("Normal")
                elif sample.endswith("T"):
                    sample_type.append("Tumor")
                else:
                    sample_type.append("Unknown")

            type = "Unknown"
            for fusion in self.__fusion_dict[key]:
                tmp = []
                for attr in fusion.zone4_attrs:
                    tmp.append(fusion.__dict__[attr])
                gene_order.append(" ".join(tmp))
                if "SameGene" in tmp:
                    type = "SameGene"
                    n_sg += 1
                elif "ReadThrough" in tmp:
                    type = "ReadThrough" if type == "Unknown" else type
                    n_rt += 1
                elif "GeneFusion" in tmp:
                    type = "GeneFusion" if type == "Unknown" else type
                    n_gf += 1
                elif "TruncatedCoding" in tmp:
                    type = "TruncatedCoding" if type == "Unknown" else type
                    n_tc +=1
                elif "TruncatedNoncoding" in tmp:
                    type = "TruncatedNoncoding" if type == "Unknown" else type
                    n_tn += 1
                elif "NoDriverGene" in tmp:
                    type = "NoDriverGene" if type == "Unknown" else type
                    n_ns += 1
                else:
                    print >> sys.stderr, "Fusions without category:", fusion.tostring()
            if type != "Unknown":
                print "Fusion", fusion.t_gene1, fusion.t_gene2, ",".join(list(set(self.__fusion_samples_dict[key]))), type, ",".join(sorted(list(set(sample_type)))),
                print "|".join(list(set(gene_order)))
                print "\tSameGene:", n_sg
                print "\tReadThrough:", n_rt
                print "\tGeneFusion:", n_gf
                print "\tTruncatedCoding:", n_tc
                print "\tTruncatedNoncoding:", n_tn
                print "\tNonSense:", n_ns


class CffFusion():
    # chr1 pos1 strand1 chr2 pos2 strand2 library sample_name sample_type disease tool split_cnt span_cnt tool_gene1 tool_gene_area1 tool_gene2 tool_gene_area2 fwd_gene_order bwd_gene_order
    def __init__(self, cff_line):
        tmp = cff_line.split()
        self.line = cff_line
        # Breadkpoint Zone
        self.chr1 = tmp[0]
        try: self.pos1 = int(tmp[1])
        except: raise ValueError("CFF Column pos1 value " + tmp[1] + " is not a valid integer\nInvalid entry: " + cff_line)
        self.strand1 = tmp[2]
        self.chr2 = tmp[3]
        try: self.pos2 = int(tmp[4])
        except: raise ValueError("CFF Column pos1 value " + tmp[4] + " is not a valid integer\nInvalid entry: " + cff_line)
        self.strand2 = tmp[5]
        # Sample info Zone
        self.library = tmp[6] # DNA/RNA
        self.sample_name = tmp[7]
        self.sample_type = tmp[8] # Tumor/Normal
        # Check to make sure sample_type is in ["Tumor", "Normal", "NA"]
        if self.sample_type not in ["Tumor","Normal", "NA"]:
          raise ValueError("sample_type value '" + tmp[8] + "' must be Tumor or Normal\nInvalid entry: " + cff_line) 
        self.disease = tmp[9]
        # Software Zone
        self.tool = tmp[10]
        try:
          self.split_cnt = int(tmp[11])
          self.span_cnt = int(tmp[12])
        except: 
          raise ValueError("split_cnt and span_cnt must be integers. Use -1 for null values\nInvalid entry: " + cff_line)
        self.t_gene1 = tmp[13] # gene reported by tool
        self.t_area1 = tmp[14] # exon/utr/intron
        self.t_gene2 = tmp[15]
        self.t_area2 = tmp[16]
        #FOR USE IN annotate_called_fusion_file.py
        self.left = []
        self.right = []
        # Re-annotation Zone
        # ReadThrough     DTX2    cds     DTX2P1-UPK3BP1-PMS2P11  utr3    True    TrueTrue     True    5.5     1       474827  1       F00000001       CCTCCCGCAGGGCCCTGAGCACCCCAATCCCGGAAAGCCGTTCACTGCCAGAGGGTTTCCCCGCCAGTGCTACCTTCCAGACAACGCCCAGGGCCGCAAG    CCTCCAGGGGCTTCCAGAACCCGGAGACACTGGCTGACATTCCGGCCTCCCCACAGCTGCTGACCGATGGCCACTACATGACGCTGCCCGTGTCTCCGGA
        if len(tmp) >= 33 :
            self.category = tmp[17]
            self.reann_gene1 = tmp[18]
            self.reann_type1 = tmp[19]
            self.reann_gene2 = tmp[20]
            self.reann_type2 = tmp[21]
            self.gene1_on_bndry = tmp[22] # gene1 breakpoint on boundary
            self.gene1_close_to_bndry = tmp[23] # gene1 breakpoint within 10bp of a boundary
            self.gene2_on_bndry = tmp[24] # gene2 on boundary
            self.gene2_close_to_bndry = tmp[25] # gene2 close to boundary
            self.score = tmp[26]    # score for bpann1 + bpann2
            self.coding_id_distance = tmp[27]   # difference between two coding fusion genes, if not coding the value is -1
            self.gene_interval_distance = tmp[28] # distance between two fusion gene intervals
            self.dnasupp = tmp[29]
            self.fusion_id = tmp[30]
            self.seq1 = tmp[31]
            self.seq2 = tmp[32]
            self.is_inframe = False
            #self.splice_site1 = "NA"
            #self.splice_site2 = "NA"
            self.closest_exon1 = tmp[34] 
            self.closest_exon2 = tmp[35]
            self.captured_reads = int(tmp[36]) 
            self.transcript1 = tmp[37]
            self.transcript2 = tmp[38]
        else:
            self.category = "NA"    # category
            self.reann_gene1 = "NA"
            self.reann_type1 = "NA"
            self.reann_gene2 = "NA"
            self.reann_type2 = "NA"
            self.gene1_on_bndry = "NA" # gene1 breakpoint on boundary
            self.gene1_close_to_bndry = "NA" # gene1 breakpoint within 10bp of a boundary
            self.gene2_on_bndry = "NA" # gene2 on boundary
            self.gene2_close_to_bndry = "NA" # gene2 close to boundary
            self.score = 0  # score for bpann1 + bpann2
            self.coding_id_distance = -1    # difference between two coding fusion genes, if not coding the value is -1
            self.gene_interval_distance = 0 # distance between two fusion gene intervals
            self.fusion_id = "NA"
            self.seq1 = "NA"
            self.seq2 = "NA"
            self.is_inframe = False
            #self.splice_site1 = "NA"
            #self.splice_site2 = "NA"
            self.closest_exon1 = "NA"
            self.closest_exon2 = "NA"
            self.transcript1 = "NA"
            self.transcript2 = "NA"
            self.captured_reads = -1
            if len(tmp) == 30:
                self.dnasupp = tmp[29]
            else:
                self.dnasupp = -9 # not available
        self.bpann1 = GeneBed("")   # breakpoint1 annotation (GeneBed)
        self.bpann2 = GeneBed("")   # breakpoint2 annotation (GeneBed)  
        self.score1 = 0 # score for bpann1 + bpann2
        self.score2 = 0 # score for bpann1 + bpann2
        
        self.boundary_info = ""     
        # same all attrs in a list, for printing    
        self.zone1_attrs = ["chr1", "pos1", "strand1", "chr2", "pos2", "strand2"]
        self.zone2_attrs = ["library", "sample_name", "sample_type", "disease"]
        self.zone3_attrs = ["tool", "split_cnt", "span_cnt", "t_gene1", "t_area1", "t_gene2", "t_area2"]
        #self.zone4_attrs = ["category", "reann_gene1", "reann_type1", "reann_gene2", "reann_type2", "gene1_on_bdry", "gene1_close_to_bndry", "gene2_on_bdry", "gene2_close_to_bndry", "score", "coding_id_distance", "gene_interval_distance", "dna_support", "fusion_id", "seq1", "seq2", "is_inframe", "splice_site1", "splice_site2", "captured_reads"]  
        self.zone4_attrs = ["category", "reann_gene1", "reann_type1", "reann_gene2", "reann_type2", "gene1_on_bdry", "gene1_close_to_bndry", "gene2_on_bdry", "gene2_close_to_bndry", "score", "coding_id_distance", "gene_interval_distance", "dna_support", "fusion_id", "seq1", "seq2", "is_inframe", "closest_exon1", "closest_exon2", "captured_reads","transcript1","transcript2"]
        #self.zone4_attrs = ["reann_gene_order1", "reann_gene_type1", "reann_gene_index1", "reann_category1", "reann_gene_order2", "reann_gene_type2", "reann_gene_index2", "reann_category2"]
        # format chr
        if not self.chr1.startswith("chr"):
            self.chr1 = "chr" + self.chr1
        if not self.chr2.startswith("chr"):
            self.chr2 = "chr" + self.chr2
        # check fields
        if not self.check_cff():
            sys.exit(1)

    def get_gene_names_from_gene_order(self):
        # gene name list
        g1 = []
        g2 = []
        g3 = []
        g4 = []
        if self.reann_gene_order1 != "NA":
            tmp = self.reann_gene_order1.split(">>")
            tmp2 = tmp[0].split(";")
            for t2 in tmp2:
                tmp3 = t2.split("_")
                g1.append(tmp3[0])
            tmp2 = tmp[1].split(";")
            for t2 in tmp2:
                tmp3 = t2.split("_")
                g2.append(tmp3[0])
            
            
        if self.reann_gene_order2 != "NA":
            tmp = self.reann_gene_order2.split(">>")
            tmp2 = tmp[0].split(";")
            for t2 in tmp2:
                tmp3 = t2.split("_")
                g3.append(tmp3[0])
            tmp2 = tmp[1].split(";")
            for t2 in tmp2:
                tmp3 = t2.split("_")
                g4.append(tmp3[0])
        
        return g1, g2, g3, g4       

    def __test_prefered_gene(self, genes, gene_name_list):
        prefered_gene = ""
        for gene in genes:
            if gene.gene_name in gene_name_list:
                if not prefered_gene:
                    prefered_gene = gene
                else:
                    if gene.type == "cds":
                        prefered_gene = gene
                        break
                    elif gene.type == "intron" and prefered_gene != "cds":
                        prefered_gene = gene
                        
        return prefered_gene        
            
    # among all the re-annotated genes, return a prefered one which should be: 1. breakpoint in/on exon. 2. if not, has a exon closest to breakpoint. 3. the longest isoform.
    # Not finished, needs to consider
    def get_perfered_ann(self, gene_ann):
        print >> sys.stderr,  "Not finished."
        return
        genes1 = gene_ann.map_pos_to_genes(self.chr1, self.pos1)
        genes2 = gene_ann.map_pos_to_genes(self.chr2, self.pos2)
        # get gene name list from re-annotated gene orders
        g1, g2, g3, g4 = self.get_gene_names_from_gene_order()
        
        prefered_gene1 = self.__test_prefered_gene(genes1, g1)
        prefered_gene2 = self.__test_prefered_gene(genes2, g2)
        # try the opposite, because breakpoints order are not always consistant with gene orders
        if not (prefered_gene1 and prefered_gene2):
            prefered_gene1 = self.__test_prefered_gene(genes1, g2)
            prefered_gene2 = self.__test_prefered_gene(genes2, g1)
            
        prefered_gene3 = self.__test_prefered_gene(genes1, g3)
        prefered_gene4 = self.__test_prefered_gene(genes2, g4)
            
        if not (prefered_gene3 and prefered_gene4):
            prefered_gene3 = self.__test_prefered_gene(genes1, g4)
            prefered_gene4 = self.__test_prefered_gene(genes2, g3)
        
        print self.tostring()
        if prefered_gene1 and prefered_gene2:
            print "\t" + prefered_gene1.tostring()
            print "\t" + prefered_gene2.tostring()
        elif prefered_gene3 and prefered_gene4:     
            print "\t" + prefered_gene3.tostring()
            print "\t" + prefered_gene4.tostring()
        else:
            print "No prefered genes."
        
    # compare fusion breakpoints and return in an order of smaller to  bigger
    def get_ordered_breakpoints(self):
        # MODIFIED, SINCE PREVIOUSLY "self.chr1 < self.chr2" was ""self.chr1 < self.chr1", which is clearly a mistake
        if self.chr1 < self.chr2:
            small_bp = (self.chr1, self.pos1, self.strand1)
            big_bp = (self.chr2, self.pos2, self.strand2)
        elif self.chr1 == self.chr2 and self.pos1 < self.pos2:
            small_bp = (self.chr1, self.pos1, self.strand1)
            big_bp = (self.chr2, self.pos2, self.strand2)
        else:
            small_bp = (self.chr2, self.pos2, self.strand2)
            big_bp = (self.chr1, self.pos1, self.strand1)
        return small_bp, big_bp 
                
    # after reannotation, get a list of upstream genes and a list of downstream genes, can be used to compare fusions from different tools on gene leverl
    def get_reannotated_genes(self, id):
        up_genes = []
        down_genes = []
        if id == 1:
            gene_order = self.reann_gene_order1
        elif id == 2:
            gene_order = self.reann_gene_order2
        else:
            print >> sys.stderr, "Strand has to be 1 or 2", id, "provided."
            sys.exit(1)
            
        #LINC00875_utr5,intron>>NBPF9_intron;LOC100288142_intron;NBPF8_intron    NoncodingGene>>CodingGene,CodingGene,CodingGene >>566_f,525_r,565_f     TruncatedNoncoding
        if gene_order != "NA":
            tmp = gene_order.split(">>")
            for g in tmp[0].split(";"):
                gname = g.split("_")[0]
                if g: up_genes.append(gname)
            for g in tmp[1].split(";"):
                gname = g.split("_")[0]
                if g: down_genes.append(gname)
        return set(up_genes), set(down_genes)

    def check_cff(self):
        if self.library not in ["NA", "DNA", "RNA"]:
            print >> sys.stderr, "Unknown library type:", self.library
            return False
        if self.sample_type not in ["NA", "Tumor", "Normal"]:
            print >> sys.stderr, "Unknown sample type:", self.sample_type
            return False
        return True

    def tostring(self):
        value = []
        for attr in self.zone1_attrs + self.zone2_attrs + self.zone3_attrs:
            if not attr in self.__dict__:
                print >> sys.stderr, "Attribute name error:", attr
                sys.exit(1)
            else:
                value.append(self.__dict__[attr])
        self.boundary_info = "\t".join(map(str, [self.gene1_on_bndry, self.gene1_close_to_bndry, self.gene2_on_bndry, self.gene2_close_to_bndry]))
        if self.fusion_id != "NA":
            return "\t".join(map(lambda x:str(x), value)) + "\t" + self.category + "\t" + self.reann_gene1 + "\t" + self.reann_type1 + "\t" + self.reann_gene2 + "\t" + self.reann_type2 + "\t" + self.boundary_info + "\t" + str(self.score) + "\t" + str(self.coding_id_distance) + "\t" + str(self.gene_interval_distance) + "\t" + str(self.dnasupp) + "\t" + self.fusion_id + "\t" + self.seq1 + "\t" + self.seq2 + "\t" + str(self.is_inframe) + "\t" + self.closest_exon1 + "\t" + self.closest_exon2 + "\t" + str(self.captured_reads) + "\t" + str(self.transcript1) + "\t"  + str(self.transcript2)
        else:
            return "\t".join(map(lambda x:str(x), value)) 
    
    #return T/F if bp on a boundary or within 5 bp of boundary designated in gene bed file
    def __check_boundary(self, bpann, order): # bpann is GeneBed object, order = head/tail
        # set on boundary info, used boundaries according to head/tail gene and their strand
        is_on_boundary = False
        close_to_boundary = False
        t = 5
    
        if bpann.start <= self.pos1 <= bpann.end:
            pos = self.pos1
        elif bpann.start <= self.pos2 <= bpann.end:
            pos = self.pos2
        else:
            return is_on_boundary, close_to_boundary
            
        
        if order == "head":
            if bpann.strand == "f": # + strand head gene
                if bpann.end == pos:
                    is_on_boundary = True
                    close_to_boundary = True
                elif abs(bpann.end-pos) < t:
                    close_to_boundary = True
            elif bpann.strand == "r":
                if bpann.start == pos:
                    is_on_boundary = True
                    close_to_boundary = True
                elif abs(bpann.start-pos) < t:
                    close_to_boundary = True
            else:
                print >> "Unknkown strand:", bpann.strand
                sys.exit(1)
        else:
            if bpann.strand == "f": # + strand head gene
                if bpann.start == pos:
                    is_on_boundary = True
                    close_to_boundary = True
                elif abs(bpann.start-pos) < t:
                    close_to_boundary = True
            elif bpann.strand == "r":
                if bpann.end == pos:
                    is_on_boundary = True
                    close_to_boundary = True
                elif abs(bpann.end-pos) < t:
                    close_to_boundary = True
            else:
                print >> "Unknkown strand:", bpann.strand
                sys.exit(1)
        return is_on_boundary, close_to_boundary

    # assign a score for every potential fusion gene pair;order is head/tail gene       
    def __cal_score(self, bpann, order):
        # Scores to choose best annotation
        score_exon_bnd = 3
        score_exon_bnd_close = 2.9
        score_utr_bnd = 2.5
        score_utr_bnd_close = 2.4
        score_in_exon = 2
        score_in_utr = 1.5
        score_in_intron = 1
        score_intergenic = 0.5
        # If gene_name not in gene_bed file, assign 0.5 as score
        if bpann.gene_name == "NA":
            score = score_intergenic
            return score
        # If gene name is in gene bed file, check boundary
        # Score according to key above
        # Exon > UTR > intron
        # On_boundary > close_to_boundary > inside_region
        is_on_boundary, close_to_boundary = self.__check_boundary(bpann, order)
        score = 0
        if bpann.type == "cds":
            if is_on_boundary:
                score = score_exon_bnd
            elif close_to_boundary:
                score = score_exon_bnd_close
            else:
                score = score_in_exon
        elif bpann.type == "utr3" or bpann.type == "utr5":
            if is_on_boundary:
                score = score_utr_bnd
            elif close_to_boundary:
                score = score_utr_bnd_close
            else:
                score = score_in_utr
        elif bpann.type == "intron":
            score = score_in_intron
        else:
            print >> sys.stderr, "Unknown type:", bpann.type
            sys.exit(1)
        return score, is_on_boundary, close_to_boundary 
                        
    # according to fusion strand (defuse style, strands are supporting pairs') return all possible gene fusions;
    # depending on gene1 and gene2 (a,b,c,d), may have to switch pos order
    # Switching should theoretically ONLY occur when a gene isnt in the gene bed file
    # Metafusion wont know how to annotate the gene if the loc or name is not in the gene_bed file
    # Arriba, starfusion and Fusioncatcher all output the genes in the correct order
    def __check_gene_pairs(self, genes1, genes2, gene_ann, switch_pos):
        # check to make sure both genes1 and genes2 have items in them
        gene_order = []
        type1 = []
        type2 = []
        id1 = []
        id2 = []
        category = ""
        t = 5
                
        for gene_name in set(genes1.keys()):
            if gene_ann.is_coding(gene_name):
                type1.append("CodingGene")
                id1.append(str(gene_ann.get_coding_gene_idx(gene_name)))
            else:
                type1.append("NoncodingGene")
        for gene_name in set(genes2.keys()):
            if gene_ann.is_coding(gene_name):
                type2.append("CodingGene")
                id2.append(str(gene_ann.get_coding_gene_idx(gene_name)))
            else:
                type2.append("NoncodingGene")
        # Calculate score for every gene pair, choose the best
        
        # get best score for current gene combination   
        max_t1 = 0, "NA", "NA", GeneBed("")
        max_t2 = 0, "NA", "NA", GeneBed("")
        # This step will pick which gene name to assign to reann_gene
        # All Metafusion annotations are based on reann_gene
        # For each gene in the fusion, identify the highest scoring transcript/gene/region 
        # Where highest scoring regions are coding on boundary bps
        # Chance for reann gene to be different from renamed gene
        # Identify region/transcript with highest score  
        # If there is more than 1 transcript/gene/region with same score, will return the FIRST occurance of the score
        # This therefore can be somewhat random if there are a lot of transcripts/genes/regions with the same score
        for gname1 in genes1:
            for bpann1 in genes1[gname1]:
                score1, is_on_boundary1, close_to_boundary1 = self.__cal_score(bpann1, "head")  
                if score1 == max_t1[0] and gname1==self.t_gene1:
                    max_t1 = score1, is_on_boundary1, close_to_boundary1, bpann1
                if score1 > max_t1[0]:
                    max_t1 = score1, is_on_boundary1, close_to_boundary1, bpann1
        for gname2 in genes2:
            for bpann2 in genes2[gname2]:
                score2, is_on_boundary2, close_to_boundary2 = self.__cal_score(bpann2, "tail")  
                #print(score2, is_on_boundary2, close_to_boundary2)
                if score2 == max_t2[0] and gname2==self.t_gene2:
                    max_t2 = score2, is_on_boundary2, close_to_boundary2, bpann2
                elif score2 > max_t2[0]:
                    max_t2 = score2, is_on_boundary2, close_to_boundary2, bpann2

        # If max_t1 + max_t2 is greater than previous score (should be zero on initialization, but if checking strand order a/b or d/c will have a previous score)
        # Assign max_t1 and max_t1 information to fusion
        if max_t1[0] + max_t2[0] > self.score1 + self.score2:
            self.score1, self.gene1_on_bndry, self.gene1_close_to_bndry, self.bpann1 = max_t1
            self.score2, self.gene2_on_bndry, self.gene2_close_to_bndry, self.bpann2 = max_t2
            self.score = self.score1 + self.score2
            self.reann_gene1 = self.bpann1.gene_name
            self.reann_gene2 = self.bpann2.gene_name
            self.reann_type1 = self.bpann1.type
            self.reann_type2 = self.bpann2.type
            ## Checking if both genes are coding genes
            if id1 and id2:
                if len(set(id1)) != 1 or len(set(id2))!= 1:
                    print >> sys.stderr, "coding id err."
                    print >> sys.stderr, id1, id2
                    sys.exit(1)
                # How many coding genes are "between" gene1 and gene2
                # IDX based on the order of the gene_bed you put in. 
                # First coding gene in gene_bed == idx_1
                # Second coding gene in gene_bed == idx_2
                # I believe v75_gene.bed is chr1->chr21->chrX->chrY
                idx1 = int(id1[0].split("_")[0])
                idx2 = int(id2[0].split("_")[0])
                self.coding_id_distance = abs(idx1 - idx2)
            gene_interval1 = gene_ann.get_gene_interval(self.bpann1.gene_name)
            gene_interval2 = gene_ann.get_gene_interval(self.bpann2.gene_name)

            # If genes are on same chromosome, how far apart are they
            if gene_interval1 and gene_interval2:
                if gene_interval1.chr == gene_interval2.chr:
                    self.gene_interval_distance = max(gene_interval1.start, gene_interval2.start) - min(gene_interval1.end, gene_interval2.end)

            ## assign fusion to a category according to best score gene pair
            # No driver gene
            # If a gene/loc not in gene bed will occur ,as well as for biological reason 
            if not genes1:
                category = "NoDriverGene"
            # map to same gene
            elif self.reann_gene1 == self.reann_gene2:
                category = "SameGene"
            else:
                # category fusions into: read through, gene fusion, truncated coding, truncated noncoding, nonsense
                # If a gene/loc not in gene bed can also get truncated coding/truncated noncoding
                gene1_is_coding = gene_ann.is_coding(self.reann_gene1)
                
                if genes2:
                    gene2_is_coding = gene_ann.is_coding(self.reann_gene2)
                else:
                    gene2_is_coding = False
                if debug: 
                    print("gene1_is_coding:", gene1_is_coding, "gene2_is_coding", gene2_is_coding)
                    #print("gene1_is_coding:", gene1_is_coding, "gene2_is_coding", gene2_is_coding)
                    #print("self.reann_gene1", self.reann_gene1, "self.reann_gene2", self.reann_gene2)
                    #print("GeneAnnotation coding gene list", gene_ann._GeneAnnotation__coding_gene_list)#[reann_gene1])
                    #print("gene_ann", gene_ann._GeneAnnotation__gene_intervals)#[reann_gene1])
                    #print("gene 1 is coding", gene_ann._GeneAnnotation__gene_intervals[self.reann_gene1].is_coding)#[reann_gene1])
                    #print("gene 1 is coding", gene_ann._GeneAnnotation__gene_intervals[self.reann_gene1].tostring())#[reann_gene1])
                    #print("gene_ann", dir(gene_ann))#[reann_gene1])
                if gene1_is_coding and gene2_is_coding:
                    for id in id1:
                        tmp = id.split("_")
                        idx = int(tmp[0])
                        strand = tmp[1]
                        #ReadThrough: gene1 and gene2 are adjacent genes (id1 - id2 = 1) or overlapping genes (id1 = id2) but breakpoints cannot map to same gene
                        if (strand == "f" and  str(idx+1) + "_f" in id2 ) or (strand == "r" and  str(idx-1) + "_r" in id2) or (id in id2):
                            category = "ReadThrough"
                        else:
                            category = "GeneFusion"
                elif gene1_is_coding:
                    category = "TruncatedCoding"
                elif not gene1_is_coding:
                    category = "TruncatedNoncoding"
                else:
                    print >> sys.stderr, "Warning: Unknown category."
                    print >> sys.stderr, type1, type2
            self.category = category
            # Return transcript assigned
            self.transcript1 = self.bpann1.transcript_id
            self.transcript2 = self.bpann2.transcript_id
        return ""

    # based on given gene annotations re-annotate cff fusions, infer possible up/downstream genes, try to fill in strand if info missing
    def ann_gene_order(self, gene_ann):
        gene_order = []
        
        # get genes which correspond to 5' and 3' breakpoints (i.e. pos1 and pos2, respectively)
        # return a list of GeneBed objects corresponding to intersecting genes 
        matched_genes1 = gene_ann.map_pos_to_genes(self.chr1, self.pos1)
        # If t_gene1 (original name assigned by caller) is in matched_genes list, no need to reannotate
        try: 
            idx=[gene.gene_name for gene in matched_genes1].index(self.t_gene1)
            matched_genes1=[matched_genes1[idx]]
        except ValueError: pass 

        matched_genes2 = gene_ann.map_pos_to_genes(self.chr2, self.pos2)
        try: 
            idx=[gene.gene_name for gene in matched_genes2].index(self.t_gene2)
            matched_genes2=[matched_genes2[idx]]
        except ValueError: pass 


        a = {} # forward strand gene at pos1
        c = {} # backward strand gene at pos1
        b = {} # forward strand gene at pos2
        d = {} # backward strand gene at pos2

        #Defuse has special stuff commented out for now as forte does not use defuse
        #if self.tool!="defuse":
        # if self.tool=="defuse":
        #     if self.strand2 == "-":
        #         self.strand2 = "+"
        #     elif self.strand2 == "+":
        #         self.strand2 = "-"
       
        for gene in matched_genes1:
            if gene.strand == "f":
                a.setdefault(gene.gene_name, []).append(gene)
            else:
                c.setdefault(gene.gene_name, []).append(gene)
        for gene in matched_genes2:
            if gene.strand == "f":
                b.setdefault(gene.gene_name, []).append(gene)
            else:
                d.setdefault(gene.gene_name, []).append(gene)

        # for tools do not provide defuse-style strand, regenerate strands, this is assuming that gene1 = 5 prime gene and gene2 = 3 prime gene
        # This seems to be doing something with DEfuse again
        # I do not believe our callers will enter this if statment.
        if self.strand1 == "NA" or self.strand2 == "NA":
            gene_interval1 = ""
            gene_interval2 = ""
            for sep in [",", "/"]:
                for gene_name in (self.t_gene1).split(sep):
                    if not gene_interval1 or gene_interval1.gene_name == "NA":
                        gene_interval1 = gene_ann.get_gene_interval(gene_name)
                    else:
                        break
                for gene_name in (self.t_gene2).split(sep):
                    if not gene_interval2 or gene_interval2.gene_name == "NA":
                        gene_interval2 = gene_ann.get_gene_interval(gene_name)
                    else:
                        break
                
            # Gene Strand : Defuse Strand
            # + >> + : + -
            # + >> - : + +
            # - >> + : - -
            # - >> - : - +
            if gene_interval1 and gene_interval2:
                if gene_interval1.strand == "f":
                    self.strand1 = "+"
                else:
                    self.strand1 = "-"
                if gene_interval2.strand == "f":
                    #self.strand2 = "-"
                    self.strand2 = "+"
                else:
                    #self.strand2 = "+"
                    self.strand2 = "-"
            else:
                # failed to fill in strand, return list with warnning info
                #gene_order.append("Strand_filling_failed")
                return gene_order

                
        # gene_order includes: 5' gene >> 3' gene, 5' gene type >> 3' gene type, 5' coding gene idx >> 3' coding gene inx, category
        # Attempts to ASSIGN SCORE, TRANSCRIPT, FUSION TYPE, and REANN GENE according to gene order 
        # Coding genes are more likely to be seleceted with higher scores
        if self.strand1 == "+" and self.strand2 == "-":
            gene_order = self.__check_gene_pairs(a, d, gene_ann, False)
            gene_order += self.__check_gene_pairs(b, c, gene_ann, True)
        elif self.strand1 == "+" and self.strand2 == "+":
            gene_order = self.__check_gene_pairs(a, b, gene_ann, False)
            gene_order += self.__check_gene_pairs(d, c, gene_ann, True)
        elif self.strand1 == "-" and self.strand2 == "-":
            gene_order = self.__check_gene_pairs(c, d, gene_ann, False)
            gene_order += self.__check_gene_pairs(b, a, gene_ann, True)
        elif self.strand1 == "-" and self.strand2 == "+":
            gene_order = self.__check_gene_pairs(c, b, gene_ann, False)
            gene_order += self.__check_gene_pairs(d, a, gene_ann, True)
  

    # realign breakpoints of this fusion to the left most, not finished, how to define "left" when genes are on different chrs 
    def left_aln_fusion_bp(self, refs):
        # provided reference file lacks fusion chr
        if not (self.chr1 in refs.references and self.chr2 in refs.references):
            return (-1, -1)
        rlen = 10
        
        #pysam use 0-based coordinates, cff use 1-based coordinates
        if self.strand1 == "+" and self.strand2 == "+":
            up_seq = refs.fetch(self.chr1, self.pos1-rlen, self.pos1+rlen)
            down_seq = sequtils.rc_seq(refs.fetch(self.chr2, self.pos2-rlen, self.pos2+rlen), "rc")
        elif self.strand1 == "+" and  self.strand2 == "-":
            up_seq = refs.fetch(self.chr1, self.pos1-rlen, self.pos1+rlen)
            down_seq = refs.fetch(self.chr2, self.pos2-rlen-1, self.pos2-1+rlen)
        elif self.strand1 == "-" and  self.strand2 == "+":
            down_seq = refs.fetch(self.chr1, self.pos1-1-rlen, self.pos1+rlen-1)
            up_seq = refs.fetch(self.chr2, self.pos2-rlen, self.pos2+rlen)
        elif self.strand1 == "-" and  self.strand2 == "-":
            down_seq = refs.fetch(self.chr1, self.pos1-rlen-1, self.pos1+rlen-1)
            up_seq = sequtils.rc_seq(refs.fetch(self.chr2, self.pos2-rlen-1, self.pos2+rlen-1), "rc")
        else:
            print >> sys.stderr, "Unknown strand:", self.strand1, self.strand2
            sys.exit(1)
        
        if len(up_seq) < rlen or len(down_seq) < rlen:
            print >> sys.stderr, "Warning: reference sequence cannot be fetched."
            print >> sys.stderr, self.tostring()
        
            return (-1, -1)

        i = rlen - 1
        while i >= 0:
            if up_seq[i].upper() == down_seq[i].upper():
                i -= 1
            else:
                break
        j = rlen
        while j < 2*rlen:
            if up_seq[j].upper() == down_seq[j].upper():
                j += 1
            else:
                break

        print up_seq.upper()[0:rlen], up_seq.lower()[rlen:]
        print down_seq.lower()[0:rlen], down_seq.upper()[rlen:]
        return (rlen-1-i, j-rlen)

    # aln fusion breakpoint to exon/utr boundary
    def aln_seq_to_boundary(self, seq1, splice_seq1, splice_seq2, seq2, distance_to_bnd, strand):
        '''
        print "seq1:", seq1[-10:]
        print "splice seq1:", splice_seq1
        print "splice seq2:", splice_seq2
        print "seq2:", seq2[:10]
        '''
        step = 0 
        gt_ag_step = -1
        if distance_to_bnd > 0: # move right
            for c in seq2:
                if step == len(splice_seq1):
                    break
                splice_site1 = splice_seq1[step:step+2]
                splice_site2 = splice_seq2[len(splice_seq2)-2+step:] + seq2[max(0, step-2):step]
                if splice_site1 == "GT" and splice_site2 == "AG":
                    gt_ag_step = step
                if step == distance_to_bnd:
                    break
                if c == splice_seq1[step]:
                    step += 1
                else:
                    break

        else: # move left
            for c in seq1[::-1]:
                if step == len(splice_seq2):
                    break
                splice_site1 = seq1[len(seq1)-step:len(seq1)-step+2] + splice_seq1[:max(0, 2-step)]
                splice_site2 = splice_seq2[-2-step:len(splice_seq2)-step]
                if splice_site1 == "GT" and splice_site2 == "AG":
                    gt_ag_step = step
                if step == abs(distance_to_bnd):
                    break
                if c == splice_seq2[-1-step]:
                    step += 1
                else:
                    break
        if step == abs(distance_to_bnd): # shift to boundary
            self.splice_site1 = splice_site1
            self.splice_site2 = splice_site2
            if strand == "f":
                self.pos1 += distance_to_bnd
                self.pos2 += distance_to_bnd
            else:
            
                self.pos1 -= distance_to_bnd
                self.pos2 -= distance_to_bnd
        elif gt_ag_step > 0: #shift to GT-AG breakpoint
            step = gt_ag_step
            self.splice_site1 = splice_seq1[step:step+2] 
            self.splice_site2 = splice_seq2[len(splice_seq2)-2+step:] + seq2[max(0, step-2):step]
            if distance_to_bnd < 0: gt_ag_step = -gt_ag_step
            if strand == "f":
                self.pos1 += gt_ag_step
                self.pos2 += gt_ag_step
            else:
            
                self.pos1 -= gt_ag_step
                self.pos2 -= gt_ag_step
        else: # do not shift
            step = 0
            self.splice_site1 = splice_seq1[step:step+2] 
            self.splice_site2 = splice_seq2[len(splice_seq2)-2+step:] + seq2[max(0, step-2):step]
            
    # get coding sequence up/downstream of fusion breakpoint. trans_ann_list is a list of GeneBed for the same transcript. ref is an opened pysam fasta
    def get_transcript_seq(self, trans_ann_list, bp, ref, gene_order):
        if not trans_ann_list:
            print "empty trans_ann_list"
            return "", "", -1, -100, "", ""
        seqs = []
        is_fw = (trans_ann_list[0].strand == "f")
        is_head = (gene_order == "head")
        splice_site = "NA"
        trans_id = trans_ann_list[0].transcript_id

        bp_in_coding = False # breakpoint of head gene in coding exon

        #get last exon idx for current transcript
        exon_ann_list = filter(lambda bpann:bpann.type=="cds", trans_ann_list)
        if exon_ann_list:
            if is_fw:
                last_exon_idx = exon_ann_list[-1].idx
            else:
                last_exon_idx = exon_ann_list[0].idx
        else:
            last_exon_idx = -1
        
        last_exon_end = -1 # for tail seq, get the position of last basepair of last exon in seq2
        distance_to_bnd = -100 #distance to exon/utr boundary, if breakpoint in exon/utr it equals to bnd - bp, thus positive, otherwise negative
        
        if is_head:
            if is_fw: # head gene on forward strand
                for ann in trans_ann_list:
                    # breakpoint in current bed annotation
                    if ann.start <= bp <= ann.end:
                        seq = ref.fetch(format_chr(ann.chr, ref), ann.start-1, bp)
                        seqs.append(seq) 
                        splice_site = ref.fetch(format_chr(ann.chr, ref), bp, bp+20)
                        if ann.type == "cds" or len(seq) < 10: # breakpoint in cds, or close to cds
                            bp_in_coding = True
                        if ann.type == "intron":
                            distance_to_bnd = ann.start - bp - 1
                        else:
                            distance_to_bnd = ann.end - bp
                        bp_ann = ann
                        
                        break
                    elif ann.type == "cds":
                        seqs.append(ref.fetch(format_chr(ann.chr, ref), ann.start-1, ann.end))
            else:
                for ann in trans_ann_list:
                    # breakpoint in current bed annotation
                    if ann.start <= bp <= ann.end:
                        seq = ref.fetch(format_chr(ann.chr, ref), bp-1, ann.end)
                        seqs.append(seq) 
                        
                        splice_site = ref.fetch(format_chr(ann.chr, ref), bp-21, bp-1)
                        if ann.type == "cds" or len(seq) < 10: # breakpoint in cds, or close to cds
                            bp_in_coding = True
                        if ann.type == "intron":
                            distance_to_bnd = bp - ann.end - 1
                        else:
                            distance_to_bnd = bp - ann.start
                        bp_ann = ann
                    elif ann.start > bp and ann.type == "cds": # for now do not consider intonic breakpoint
                        seqs.append(ref.fetch(format_chr(ann.chr, ref), ann.start-1, ann.end))

        else:
            if is_fw: # tail gene on forward strand
                for ann in trans_ann_list:
                    if ann.type=="cds" and ann.idx == last_exon_idx:
                        last_exon_end = len("".join(seqs)) + (ann.end  - ann.start + 1)
                    # breakpoint in current bed annotation
                    if ann.start <= bp <= ann.end:
                        splice_site = ref.fetch(format_chr(ann.chr, ref), bp-21, bp-1)
                        seqs.append(ref.fetch(format_chr(ann.chr, ref), bp-1, ann.end)) 
                        if ann.type == "cds":
                            bp_in_coding = True
                        bp_ann = ann
                    elif ann.start > bp and ann.type != "intron": # for now do not consider intonic breakpoint
                        seqs.append(ref.fetch(format_chr(ann.chr, ref), ann.start-1, ann.end))
            else:
                for ann in trans_ann_list:
                    # breakpoint in current bed annotation
                    if ann.start <= bp <= ann.end:
                        seqs.append(ref.fetch(format_chr(ann.chr, ref), ann.start-1, bp)) 
                        splice_site = ref.fetch(format_chr(ann.chr, ref), bp, bp+20)
                        if ann.type == "cds":
                            bp_in_coding = True
                        bp_ann = ann
                        break
                    elif ann.type != "intron": # for now do not consider intonic breakpoint
                        seqs.append(ref.fetch(format_chr(ann.chr, ref), ann.start-1, ann.end))

        trans_seq = "".join(seqs).upper()
        splice_site = splice_site.upper()       
        if not is_fw:
            trans_seq = sequtils.rc_seq(trans_seq, "rc")
            splice_site = sequtils.rc_seq(splice_site, "rc")
        return trans_seq, splice_site, last_exon_end, distance_to_bnd, bp_ann, bp_in_coding 
            


    # for given cff fusion, test all combination of head/tail gene sequences to check their codons. Can be used to classiy readthrough more strictly by requiring the fusion seq in-frame and reach the end of tail gene        
    def check_codon(self, gene_ann, ref_file):
        head_seqs = []
        tail_seqs = []
        ref = pysam.FastaFile(ref_file) 
        # get all transcripts for head gene 
        head_gene_intervals = gene_ann.get_gene_interval(self.reann_gene1)
        tids = head_gene_intervals.transcript_ids
        if not head_gene_intervals.is_coding:
            print >> sys.stderr, "head gene is not coding"
            return
            
        if not tids:
            print >> sys.stderr, "Head gene not in annotation:", self.reann_gene1
            return

        for transcript_id in tids:
            bpann_list = gene_ann.get_transcripts_ann(transcript_id) # all GeneBed annotations of transcript_id
            if not bpann_list[0].start <= self.pos1 <= bpann_list[-1].end: #skip the transcript that doesnt include the breakpoint
                continue
            seq, splice_site, last_exon_end, distance_to_bnd, bp_ann, bp_in_coding =  self.get_transcript_seq(bpann_list, self.pos1, ref, "head")
            head_seqs.append((seq, splice_site, last_exon_end, transcript_id, distance_to_bnd, bp_ann, bp_in_coding))
        # get all transcripts for tail gene 
        tids = (gene_ann.get_gene_interval(self.reann_gene2)).transcript_ids
        if not tids:
            # no tail gene, make a fake bed annotation, extract 1k sequence downstream the breakpint
            fake_ann = GeneBed("")
            fake_ann_len = 1000
            fake_ann.chr = self.chr2
            fake_ann.type = "intergenic"
            fake_ann.idx = 0
            fake_ann.start = self.pos2 - fake_ann_len
            fake_ann.end = self.pos2 + fake_ann_len

            # according to head gene strand, infer fake_ann strand
            if self.strand1 == self.strand2: # head gene and tail gene on different strand   
                if head_gene_intervals.strand == "f":
                    fake_ann.strand = "r"
                else:
                    fake_ann.strand = "f"
            else: # head and tail gene on same strand
                fake_ann.strand = head_gene_intervals.strand
            seq, splice_site, last_exon_end, distance_to_bnd, bp_ann, bp_in_coding = self.get_transcript_seq([fake_ann], self.pos2, ref, "tail")
            tail_seqs.append((seq, splice_site, last_exon_end, fake_ann.transcript_id, distance_to_bnd, bp_ann, bp_in_coding))
        for transcript_id in tids:
            bpann_list = gene_ann.get_transcripts_ann(transcript_id) # all GeneBed annotations fo transcript_id
            if not bpann_list[0].start <= self.pos2 <= bpann_list[-1].end:
                continue
            seq, splice_site, last_exon_end, distance_to_bnd, bp_ann, bp_in_coding = self.get_transcript_seq(bpann_list, self.pos2, ref, "tail")
            tail_seqs.append((seq, splice_site, last_exon_end, transcript_id, distance_to_bnd, bp_ann, bp_in_coding))
        # search first stop codon in tail gene sequence, stop codon: TAA, TAG, TGA
        stop_codons= set(["TAA", "TAG", "TGA"])
        seq_set = set()
        # among all the combinations of head and tail sequences, 
        # 1. check if the fusion can reach any transcript's end (first stop codon is the tail transcript's stop codon); 
        # 2. if not 1, return the smallest stop condon position.

        # Format: infered seq, reach last exon, first stop codon pos, upstream transcript id, downstream gene transcript_id, distance_to_bnd, score
        infered_fusion_seq_info = ("", False, -1, "", "", -100, "", 0) 

        for seq1, splice_site1, last_exon_end1, transcript_id1, distance_to_bnd1, bp_ann1, bp_in_coding1 in head_seqs:
            if not seq1:
                continue
            shift1 = len(seq1) % 3
            for seq2, splice_site2, last_exon_end2, transcript_id2, distance_to_bnd2, bp_ann2, bp_in_coding2 in tail_seqs:
                i = 0
                shift_seq2 = seq2[(3-shift1)%3:]
                # convert shift_seq2 to triplets
                shift_seq2_triplets = re.findall(r'.{3}', shift_seq2)
                # get uniq triplets set
                shift_seq2_codon_set = set(shift_seq2_triplets)
                # compare to stop codons
                shift_seq2_stop_conons =  stop_codons & shift_seq2_codon_set
                # stop codon found, get the first one
                if shift_seq2_stop_conons:
                    first_stop_idx = len(shift_seq2_triplets)
                    first_stop_codon = "NA"
                    for codon in shift_seq2_stop_conons:
                        stop_idx = shift_seq2_triplets.index(codon) 
                        first_stop_idx = stop_idx
                        first_stop_codon = codon
                    tail_gene_seq = seq2[:(3-shift1)%3 + (first_stop_idx + 1)*3]    
                    update = False
                    score = 0 # current 
                    is_inframe = False
                    if bp_in_coding1 and bp_in_coding2 and (len(seq1) + len(seq2)) % 3 == 0: # head transcript is coding and inframe
                        is_inframe = True
                        score += 10000

                    if abs(infered_fusion_seq_info[5]) > abs(distance_to_bnd1): # current combination closer to annotation boundary
                        score += 1000
                    if bp_ann1.type == "cds": # head transcript invole exon
                        score += 500
                    if bp_ann2.type == "cds": 
                        score += 300
                    if bp_ann1.type == "utr5" or  bp_ann1.type == "utr3": # head transcript involve utr
                        score += 200
                    if bp_ann2.type == "utr5" or bp_ann2.type == "utr3": 
                        score += 100
                    if gene_ann.transcript_is_coding(transcript_id2): # current transcript has cds
                        score += 10
                    if first_stop_idx <= int(infered_fusion_seq_info[2]): # saved combination has a smaller stop codon position
                        score += 1
                    update = True if score > int(infered_fusion_seq_info[7]) else False
                    if update:
                
                        # infered fusion sequence, ends at first stop codon
                        infered_fusion_seq = seq1 + "_" + splice_site1 + "__" + splice_site2 + "_" + seq2[:(3-shift1)%3 + (first_stop_idx + 1)*3]
                        infered_fusion_seq_info = (infered_fusion_seq, is_inframe, first_stop_idx, transcript_id1, transcript_id2, distance_to_bnd1, bp_ann1, score)
                else:
                    pass
                    
        is_inframe = infered_fusion_seq_info[1]
        seqs = infered_fusion_seq_info[0].split("_")
        distance_to_bnd = infered_fusion_seq_info[5]
        bp_ann = infered_fusion_seq_info[6]
        if len(seqs) == 5:
            self.aln_seq_to_boundary(seqs[0], seqs[1], seqs[3], seqs[4], distance_to_bnd, bp_ann.strand)
            self.is_inframe = is_inframe
        ref.close()
        
        
class GeneInterval():
    """
    Instance variables: gene_name, chr, strand, start, end, is_coding, is_contradictory, transcript_ids
    """
    def __init__(self, bed_ann_list):
        self.__load(bed_ann_list)

    def check_bed_ann_list(self, bed_ann_list):
        flag = True
        if not bed_ann_list:
            flag = False
        else:
            if len(set([(a.gene_name, a.chr, a.strand) for a in bed_ann_list])) > 1:
                flag = False
            else:
                max_start = max([a.start for a in bed_ann_list])
                min_end = min([a.end for a in bed_ann_list])
                # Gene has two annotations more than 1Mb away from each other
                if max_start - min_end > 5000000:
                    flag = False
            if not flag:
                print >> sys.stderr, "Warning: Input gene annotations include multiple chr, strand, or regions (5Mb away). Skipping current gene annotation."
                print >> sys.stderr, set([(a.gene_name, a.chr, a.strand) for a in bed_ann_list])

        return flag

    def __load(self, bed_ann_list):
        if not self.check_bed_ann_list(bed_ann_list):
            self.gene_name = "NA"
            self.chr = "NA"
            self.start = -1
            self.end = -1
            self.strand = "NA"
            self.is_coding = False
            self.is_contradictory = True # contradictory gene annotation, will not be used
            self.transcript_ids = []
        else:
            if debug: print([gene_bed.tostring() for gene_bed in bed_ann_list if gene_bed.type == "cds"])
            self.gene_name = bed_ann_list[0].gene_name
            self.chr = bed_ann_list[0].chr
            self.strand= bed_ann_list[0].strand
            self.start = min([a.start for a in bed_ann_list])
            self.end = max([a.end for a in bed_ann_list])
            self.is_coding = True if "cds" in [a.type for a in bed_ann_list] else False
            self.is_contradictory = False
            self.transcript_ids = list(set([genebed.transcript_id for genebed in bed_ann_list]))

    def tostring(self): 
        return "\t".join([self.gene_name, self.chr, self.strand, str(self.start), str(self.end), str(self.is_coding), str(self.is_contradictory)])

    def overlap(self, interval2):
        if self.chr == interval2.chr and min(self.end, interval2.end) - max(self.start, interval2.start) > 0:
            return True
        else:
            return False
    def merge(self, interval2):
        if self.chr != interval2.chr or self.strand != interval2.strand:
            print >> sys.stderr, "Warning: intervals are on different chr/strand."
            print >> sys.stderr, self.gene_name, interval2.gene_name
            return self
        else:
            new_interval = copy.copy(interval2)
            new_interval.gene_name = "Merged_" + self.gene_name + "_" + interval2.gene_name
            new_interval.chr =  self.chr
            new_interval.start =  min(self.start, interval2.start)
            new_interval.end =  max(self.end, interval2.end)
            new_interval.strand = self.strand
            new_interval.transcript_ids = list(set(self.transcript_ids + interval2.transcript_ids))
            # for merged gene intervals is_coding has no sense
            new_interval.is_coding = False
            
            return new_interval         

class GeneBed():
    """
    GeneBed class is used to represent a line of gene_ann_bed file.
 
    :param gene_ann_bed: an ensemble gene file containing: 
	    chr	        start	        end	        transcript_id	type	idx	strand	gene_name	gene_id	
            chr1	24687531	24696163	ENST00000003583	intron	1	r	STPG1	        ENSG00000001460

    :instance variables: chr, start, end, transcript_id, type, idx, strand, gene_name, gene_id, is_from_contradictory_gene 

    """
    def __init__(self, bed_line):
        if bed_line:
            tmp = bed_line.split()
            self.chr = tmp[0]
            self.start = int(tmp[1]) + 1    # bed file use 0-based coordinate
            self.start_orig = int(tmp[1]) 
            self.end = int(tmp[2])      # start and end are first and last base of each segment
            self.transcript_id = tmp[3]
            self.type = tmp[4] # utr/intron/cds
            self.idx = int(tmp[5])
            self.strand = tmp[6] # r/f
            self.gene_name = tmp[7]
            self.gene_id = tmp[8]
            # the following two attr are only used for breakpoint annotation, tell whether current breakpoint on/close to boundary, need to reset everytime before used
            #self.is_on_boundary = False
            #self.close_to_boundary = False
            # in gene intervals, this gene includes annotations from different chr, strand or >5Mb distance
            self.is_from_contradictory_gene = False
        else:
            self.chr = "NA"
            self.start = 0  # bed file use 0-based coordinate
            self.end = 0        # start and end are first and last base of each segment
            self.transcript_id = "NA"
            self.type = "NA" # utr/intron/cds
            self.idx = -1
            self.strand = "NA"
            self.gene_name = "NA"
            self.gene_id = "NA"
            self.is_from_contradictory_gene = False

    def overlap(self, genebed2):
        if self.chr == genebed2.chr and min(self.end, genebed2.end) - max(self.start, genebed2.start) > 0:
            return True
        else:
            return False

    def merge(self, genebed2):
        if self.chr != genebed2.chr or self.strand != genebed2.strand:
            print >> sys.stderr, "Warning: intervals are on different chr/strand."
            print >> sys.stderr, self.gene_name, genebed2.gene_name
            return self
        else:
            new_bed = copy.copy(genebed2)
            new_bed.transcript_id = "Merged"
            new_bed.start =  min(self.start, genebed2.start)
            new_bed.end =  max(self.end, genebed2.end)
            
            return new_bed

    def tostring(self):
        attrs = [self.chr, str(self.start_orig), str(self.end), self.transcript_id, self.type, str(self.idx), self.strand, self.gene_name, self.gene_id]
        return "\t".join(attrs)

#Load bed format gene annotation, current support knowngene.bed's format, map given genomic loactions to genens, return matched gene list
class GeneAnnotation():
    __gene_starts = {}
    __genes = {}
    __max_diff = 1000000
    __gene_intervals = {}
    __transcripts_ann = {}
    __gene_name_id_map = {}
    __gene_name_idx_map = {}
    __coding_gene_list = []
    # gene_ann_bed
    # chr1    24683494  24685032    ENST00000003583 utr3    0       r       STPG1   ENSG00000001460
    def __init__(self, gene_ann_bed):
        """
 
        """
        if gene_ann_bed != "":
            # load gene_ann_bed file into 
            self.load_gene_bed(gene_ann_bed)
            # create dictionary of {gene_name :  GeneInterval } 
            self.load_gene_intervals(gene_ann_bed)
            # builds coding gene list, and also indexes all coding gene intervals
            self.build_coding_gene_list()
            '''
            print >> sys.stderr, "mem of gene_starts:", sys.getsizeof(self.__gene_starts)
            print >> sys.stderr, "length:", len(self.__gene_starts.keys())
            print >> sys.stderr, "unit size:", sys.getsizeof(self.__gene_starts[self.__gene_starts.keys()[0]])
                
            print >> sys.stderr, "mem of genes:", sys.getsizeof(self.__genes)
            print >> sys.stderr, "length:", len(self.__genes.keys())
            print >> sys.stderr, "unit size:", sys.getsizeof(self.__genes[self.__genes.keys()[0]])

            print >> sys.stderr, "mem of gene_intervals:", sys.getsizeof(self.__gene_intervals)
            print >> sys.stderr, "length:", len(self.__gene_intervals.keys())
            print >> sys.stderr, "unit size:", sys.getsizeof(self.__gene_intervals[self.__gene_intervals.keys()[0]].start)
            '''
            self.load_transcripts_ann(gene_ann_bed)

    def load_gene_bed(self, gene_ann_bed):
        """
        : self.__genes : A dictionary of { chr : [list of GeneBed objects] }
                         For each chromosome, lists all genomic features on that chromosome in the form of GeneBed objects
        """
        start_time = time.time()
        n = 0
        for line in open(gene_ann_bed, "r"):
            ann = GeneBed(line)
            # filter out annotation not needed, e.g. gene annotation on chr_hap
            if self.filter_gene_ann(ann):
                continue
            self.__genes.setdefault(ann.chr, []).append(ann)
            # gene id and name map, this dict will take ~3G ram for ensgene annotation, current not loaded
            #self.__gene_name_id_map.setdefault(ann.gene_name, ann.gene_id)
            n += 1
        print >> sys.stderr, n, "annotations from", gene_ann_bed, "loaded."


        # sort all ann of each chr by start pos
        for chr in self.__genes:
            self.__genes[chr] = sorted(self.__genes[chr], key=lambda d:int(d.start))

            # save start pos in self.__gene_starts
            self.__gene_starts[chr] = [int(d.start) for d in self.__genes[chr]]
        print >> sys.stderr, time.time() - start_time, "sec. elapsed."  

        # TEST...checking self.__genes dictionary
        #for key in self.__genes.keys():
        #    print key, self.__genes[key]
        #... TEST

    def load_gene_intervals(self, gene_ann_bed):
        """
        For each gene use minimal start and max end of all its genomic features as its interval, 
        save intervals of all genes in a dictionary self.__gene_intervals, which is a dict of
        {gene_name :  GeneInterval }
        """
        start_time = time.time()
        n1 = 0
        n2 = 0
        tmp_dict = {}
        # use gene name as key to bulid a dict for gene_ann_bed, all annotations for the same gene are saved in a list which will be used for getting intervals
        for line in open(gene_ann_bed, "r"):
            bed_ann = GeneBed(line)
            # filter out annotation not needed, e.g. gene annotation on chr_hap
            if self.filter_gene_ann(bed_ann):
                continue
            # builds dict of { gene_name : [list of gene features] }
            tmp_dict.setdefault(bed_ann.gene_name, []).append(bed_ann)

        # TEST...
       # for key in tmp_dict.keys():
       #     print key, [(geneBed.gene_id, geneBed.start, geneBed.end ) for geneBed in tmp_dict[key]]
            # CTD-2588J6.2 [('ENSG00000261029', 82203901, 82203952), ('ENSG00000261029', 82203953, 82209036), ('ENSG00000261029', 82209037, 82209408)]
        #... TEST
        if debug == 2: print("tmp_dict", [[gene for gene in tmp_dict[key]] for key in tmp_dict.keys()])
        #debug
        for gene_name in tmp_dict:
            self.__gene_intervals.setdefault(gene_name, GeneInterval(tmp_dict[gene_name]))

        # TEST...confirming that gene intervals are correctly acquired: CONFIRMED!
        #for key in self.__gene_intervals.keys():
        #    print key, self.__gene_intervals[key].start, self.__gene_intervals[key].end
        #... TEST
        
        print >> sys.stderr, "Gene intervals loaded."
        print >> sys.stderr, time.time() - start_time, "sec. elapsed." 
 
    def build_coding_gene_list(self):
        """
        Populates index dictionary self.__gene_name_idx_map, a dictionary of { gene_name : index_value } pairs
        Indices are assigned to each forward and reverse strand gene. Indices allows us to determine whether two genes are adjacent,
        since adjacent genes might form a readthrough fusion
        """
        #Populates self.__coding_gene_list with GeneInterval objects that represent coding genes 
        for gene_name in self.__gene_intervals:
            interval = self.__gene_intervals[gene_name]
            if interval.is_coding: #
                self.__coding_gene_list.append(interval)
        self.__coding_gene_list.sort(key = lambda i:(i.chr, i.start))

        i_f = 0 # idx of forward strand gene
        i_r = 0 # idx of reverse strand gene
        pre_interval_f = GeneInterval([])
        pre_interval_r = GeneInterval([])
        # merges  overlapping intervals into a single GeneInterval object 
        for interval in self.__coding_gene_list:
            # assign unique indices to "forward" strand gene intervals
            if interval.strand == "+" or interval.strand == "f":
                if interval.overlap(pre_interval_f):
                    pre_interval_f = interval.merge(pre_interval_f)
                else:
                    pre_interval_f = interval
                    i_f += 1
                self.__gene_name_idx_map.setdefault(interval.gene_name, str(i_f) + "_" + interval.strand)
            # assign unique indices to "reverse" strand gene intervals
            elif interval.strand == "-" or interval.strand == "r":
                if interval.overlap(pre_interval_r):
                    pre_interval_r = interval.merge(pre_interval_r)
                else:
                    pre_interval_r = interval
                    i_r += 1
                self.__gene_name_idx_map.setdefault(interval.gene_name, str(i_r) + "_" + interval.strand)

    def load_transcripts_ann(self, gene_ann_bed):
        """
        Populate dictionary __transcripts_ann, a dictionary of { transcript_id : [list of GeneBed objects] }.
        Allows GeneBed object to be retrieved using transcript_id
        """
        start_time = time.time()
        for line in open(gene_ann_bed, "r"):
            ann = GeneBed(line)
            key = ann.transcript_id
            self.__transcripts_ann.setdefault(key, []).append(ann)
        # TEST..
        #for key in self.__transcripts_ann:
        #    print key, [item for item in self.__transcripts_ann[key] ]
        #... TEST
        print >> sys.stderr, "Transcript annotations loaded."
        print >> sys.stderr, time.time() - start_time, "sec. elapsed."  
    
    def get_transcripts_ann(self, trans_id):
        if trans_id in self.__transcripts_ann:
            return self.__transcripts_ann[trans_id]
        else:
            return []

    def transcript_is_coding(self, trans_id):
        if trans_id in self.__transcripts_ann:
            if filter(lambda ann:ann.type=="cds" ,self.__transcripts_ann[trans_id]):
                return True
        return False    

    def get_previous_ann(self, bpann):
        bpann_list = self.get_transcripts_ann(bpann.transcript_id)
        pre_bpann_list = []
        for ann in bpann_list:
            if bpann.strand == "f":
                if ann.idx == bpann.idx - 1:
                    pre_bpann_list.append(ann)
            elif bpann.strand == "r":
                if ann.idx == bpann.idx + 1:
                    pre_bpann_list.append(ann)
        return pre_bpann_list

    def get_next_ann(self, bpann):
        bpann_list = self.get_transcripts_ann(bpann.transcript_id)
        next_bpann_list = []
        for ann in bpann_list:
            if bpann.strand == "f":
                if ann.idx == bpann.idx + 1:
                    next_bpann_list.append(ann)
            elif bpann.strand == "r":
                if ann.idx == bpann.idx - 1:
                    next_bpann_list.append(ann)
        return next_bpann_list

    def get_gene_id(self, gene):
        if gene in self.__gene_name_id_map:
            return self.__gene_name_id_map[gene]
        else:
            return "NA"


    def print_coding_gene_idx(self):
        for key in self.__gene_name_idx_map:
            print key, self.__gene_name_idx_map[key]
    # whether a gene include coding exon (cds)
    def is_coding(self, gene_name):
        return self.__gene_intervals[gene_name].is_coding 
    # where a gene has contradictory annotations
    def is_contradictory(self, gene_name):
        return self.__gene_intervals[gene_name].is_contradictory
    def get_coding_gene_idx(self, gene_name):
        return self.__gene_name_idx_map[gene_name]

    def test_intervals(self):
        return self.__gene_intervals
    def filter_gene_ann(self, ann):
        if "hap" in ann.chr or "_" in ann.chr:
            return True
        else:
            return False

    # for a given gene, return its interval if it is in the annotation, return empty interval if not.
    def get_gene_interval(self, gene_name):
        if gene_name in self.__gene_intervals:
            return self.__gene_intervals[gene_name]
        else:
            if gene_name != "NA":
                print >> sys.stderr, "Warnning: gene name", gene_name, "is not in current annotation."
            return GeneInterval([])

    # get distance between two gene intervals, return -1 if genes are rom different chrs, minus value if overlap.
    def get_gene_distance(self, gene1, gene2):
        interval1 = self.get_gene_interval(gene1)
        interval2 = self.get_gene_interval(gene2)
        if interval1 and interval2:
            chr1 = interval1[0]
            start1 = int(interval1[1])
            end1 = int(interval1[2])
            chr2 = interval2[0]
            start2 = int(interval2[1])
            end2 = int(interval2[2])
            if chr1 != chr2: # genes from different chrs
                return -1 
            else:
                return cmp_intervals(start1, end1, start2, end2)
            
    #return a list of GeneBed annotation        
    def map_pos_to_genes(self, chr, pos):
        # if pos is within 10bp window of any boundary, set close_to_boundary True
        t = 10
        matched_genes = []
        if not chr in self.__gene_starts:
            return matched_genes
        idx = bisect.bisect(self.__gene_starts[chr], pos)   
        while 0 < idx <= len(self.__gene_starts[chr]):
            #bpann is an GeneBed object
            bpann = self.__genes[chr][idx-1]
            #search within a limited region (default 1000000)
            if pos - bpann.start > self.__max_diff: 
                break

            if bpann.start <= pos <= bpann.end:
                matched_genes.append(bpann)
            idx -= 1
    
        return matched_genes
    def get_adjacent_exons(self, chr, pos):
        previous_exons = []
        next_exons = [] 
            
        if not chr in self.__gene_starts:
            return previous_exons, next_exons
        idx = bisect.bisect(self.__gene_starts[chr], pos)
        while 0 < idx <= len(self.__gene_starts[chr]):
            
            bpann = self.__genes[chr][idx-1]
            #search within a limited region (default 1000000)
            if pos - bpann.start > self.__max_diff: 
                break
            if bpann.start <= pos <= bpann.end:
                # when breakpoint map to an intron/exon annotation
                #if bpann.type == "intron" or bpann.type == "cds": 
                if True: 
                    # previous and next exons are based on coordinates, not strand. i.e. if a adjacent exon's coordinate < bp, it is a previous exon, otherwise next exon
                    for i in range(idx, min(idx+100, len(self.__gene_starts[chr]))):
                        adjacent_bpann = self.__genes[chr][i]
                        # in next 100 annotations, try to find one exon next to current intron of the same transcript
                        if bpann.transcript_id == adjacent_bpann.transcript_id and adjacent_bpann.type == "cds" and abs(adjacent_bpann.idx - bpann.idx) <= 1:

                            next_exons.append(adjacent_bpann)
                            #break
                            #if bpann.strand == "f":
                            #   next_exons.append(adjacent_bpann)
                            #   break
                            #else:
                            #   previous_exons.append(adjacent_bpann)
                            #   break
                    for i in range(idx-2, max(0, idx-100), -1):
                        adjacent_bpann = self.__genes[chr][i]
                        # in next 100 annotations, try to find one exon next to current intron of the same transcript
                        if bpann.transcript_id == adjacent_bpann.transcript_id and adjacent_bpann.type == "cds" and abs(adjacent_bpann.idx - bpann.idx) <= 1:
                            previous_exons.append(adjacent_bpann)
                            #break
                            #if bpann.strand == "r":
                            #   next_exons.append(adjacent_bpann)
                            #   break
                            #else:
                            #   previous_exons.append(adjacent_bpann)
                            #   break
            idx -= 1
        return previous_exons, next_exons   

    def get_closest_exon_lists(self, chr, pos):
        # returns exon if breakpoint (pos) falls within it, otherwise returns 2 exons that flank breakpoint. previous and next exons are based on coordinates, not strand. i.e. if a adjacent exon's coordinate < bp, it is a previous exon, otherwise next exon
        previous_exons = []
        next_exons = []
        if not chr in self.__gene_starts:
            raise Exception(" 'chr' prefix needed for chromosomes")
            return previous_exons, next_exons
        idx = bisect.bisect(self.__gene_starts[chr], pos)
        #print("idx", idx)
        while 0 < idx <= len(self.__gene_starts[chr]):
            bpann = self.__genes[chr][idx-1]
            #search within a limited region (default 1000000)
            if pos - bpann.start > self.__max_diff:
                break
            # search all annotations that breakpoint is within 
            if bpann.start <= pos <= bpann.end:
                # if current annnotation is a cds, we have found closest exon
                if bpann.type in ["cds", "utr5", "utr3"]:
                    return [bpann], [bpann]
                # search for "next" exon (i.e on right)
                for i in range(idx, min(idx+100, len(self.__gene_starts[chr]))):
                    adjacent_bpann = self.__genes[chr][i]
                    # in next 100 annotations, try to find exon or UTR adjacent to current intron of the same transcript
                    if bpann.transcript_id == adjacent_bpann.transcript_id and adjacent_bpann.type in ["cds", "utr5", "utr3"] and abs(adjacent_bpann.idx - bpann.idx) <= 1:
                        next_exons.append(adjacent_bpann)
                # search for "previous" exon (i.e on left)
                for i in range(idx-2, max(0, idx-100), -1):
                    adjacent_bpann = self.__genes[chr][i]
                    # in next 100 annotations, try to find exon adjacent to current intron of the same transcript
                    if bpann.transcript_id == adjacent_bpann.transcript_id and adjacent_bpann.type == ["cds", "utr5", "utr3"] and abs(adjacent_bpann.idx - bpann.idx) <= 1:
                        previous_exons.append(adjacent_bpann)
            idx -= 1
        return previous_exons, next_exons

    def get_adjacent_introns(self, chr, pos):
        previous_introns = []
        next_introns = []   
            
        if not chr in self.__gene_starts:
            return previous_introns, next_introns
        idx = bisect.bisect(self.__gene_starts[chr], pos)
        while 0 < idx <= len(self.__gene_starts[chr]):
            
            bpann = self.__genes[chr][idx-1]
            #search within a limited region (default 1000000)
            if pos - bpann.start > self.__max_diff: 
                break
            if bpann.start <= pos <= bpann.end:
                # when breakpoint map to an intron/exon annotation
                #if bpann.type == "intron" or bpann.type == "cds": 
                if True: 
                    # previous and next exons are based on coordinates, not strand. i.e. if a adjacent exon's coordinate < bp, it is a previous exon, otherwise next exon
                    for i in range(idx, min(idx+100, len(self.__gene_starts[chr]))):
                        adjacent_bpann = self.__genes[chr][i]
                        # in next 100 annotations, try to find one exon next to current intron of the same transcript
                        if bpann.transcript_id == adjacent_bpann.transcript_id and adjacent_bpann.type == "intron" and abs(adjacent_bpann.idx - bpann.idx) <= 1:
                            next_introns.append(adjacent_bpann)
                    for i in range(idx-2, max(0, idx-100), -1):
                        adjacent_bpann = self.__genes[chr][i]
                        # in next 100 annotations, try to find one exon next to current intron of the same transcript
                        if bpann.transcript_id == adjacent_bpann.transcript_id and adjacent_bpann.type == "intron" and abs(adjacent_bpann.idx - bpann.idx) <= 1:
                            previous_introns.append(adjacent_bpann)
            idx -= 1
        return previous_introns, next_introns   
                
def cmp_fusion_breakpoints(bp1, bp2, diff):
    chr1, pos1, strand1 = bp1
    chr2, pos2, strand2 = bp2
            
    if chr1 != chr2 or strand1 != strand2:
        return False
    elif abs(pos1 - pos2) < diff:
        return True
    else:
        return False    

#compare two genomic locations, return True if the difference is samller than parameter diff        
def cmp_breakpoints(chr1, pos1, chr2, pos2, diff):
    if chr1 != chr2:
        return False
    elif abs(pos1 - pos2) < diff:
        return True

# check whether two intervals overlap, overlap if end-start > 0
def cmp_intervals(start1, end1, start2, end2):
    start = max(start1, start2)
    end = min(end1, end2)
    return end - start


def pass_filter(read):
    if read.mapping_quality == 255: # tophat2 max mapq
        return True
    else:
        return False

# return #read per basepair for all exons(cds) of transcripts in trans_id_list
def get_exon_cov(trans_id_list, gene_ann, bam):

    # get uniquely mapped read number in each cds
    for trans_id in trans_id_list:
        for transcript in  gene_ann.get_transcripts_ann(trans_id):
            read_cnt = 0
            n = 0
            #if type == "cds":
            if 1:
                for read in bam.fetch(transcript.chr, transcript.start, transcript.end):
                    n += 1
                    # bam.fetch may return reads not in given region (bug of pysam?)
                    # also require read to pass a "filter", currently only requires mapping quality to be maximum, 50 for tophat2, 60 for bwa, 255 for star
                    if read.reference_start in range(transcript.start, transcript.end) and pass_filter(read):
                        read_cnt += 1
                    else:
                        pass
                        #print read.reference_start, read.reference_start in [start, end]
            
                print transcript.tostring, read_cnt, float(read_cnt)/(transcript.end - transcript.start)

def get_fusion_exon_cov(gene_fusion, gene_ann, bam):
    # map fusion breakpoints to genes
    matched_genes1 = gene_ann.map_pos_to_genes(gene_fusion.chr1, gene_fusion.pos1)  
    matched_genes2 = gene_ann.map_pos_to_genes(gene_fusion.chr2, gene_fusion.pos2)  
    
    # get all transcripts of genes
    trans_id_list1 = list(set(g[3] for g in matched_genes1))
    trans_id_list2 = list(set(g[3] for g in matched_genes2))
    trans_id_list = trans_id_list1 + trans_id_list2
    get_exon_cov(trans_id_list, gene_ann, bam)

#build fusion reference, refs is pysam.FastaFile, seg_len is the length of extraced ref sequence
def output_fusion_fasta(fusion, refs, seg_len):
    ## defuse strands have different meanings, need to convert to the following pattern before uesd
    # ++ : ->| |->
    # +- : ->| <-|
    # -+ : |<- |->
    # -- : |<- <-|
    ## defuse:
    # ++ : ->| <-|
    # +- : ->| |->
    # -+ : |<- <-|
    # -- : |<- |->
    
    # pysam.FastaFile uses 0-based coordinates, cff fusion uses 1-based
    # strand used here is differnt from defuse
    if strand1 == "+": ## ->|
        win_start1 = fusion.pos1 - seg_len
        win_end1 = fusion.pos1
        seq1 = refs.fetch(fusion.chr1, win_start1, win_end1)
    else: ## |<-
        win_start1 = fusion.pos1 - 1
        win_end1 = fusion.pos1 + seg_len - 1
        seq1 = rc_seq(refs.fetch(fusion.chr1, win_start1, win_end1))

    if strand2 == "+": ## |->
        win_start2 = fusion.pos2 - 1
        win_end2 = fusion.pos2 + seg_len - 1
        seq2 = refs.fetch(fusion.chr2, win_start2, win_end2)
    else: ## <-|
        win_start2 = fusion.pos2 - seg_len
        win_end2 = fusion.pos2
        seq2 = rc_seq(refs.fetch(fusion.chr2, win_start2, win_end2))
    seq = seq1 + seq2

# build transcript sequence based on bed format annotation file
# chr1    24683494  24685032    ENST00000003583 utr3    0       r       STPG1   ENSG00000001460 
def build_transcript_seq(gene_ann, ref):
    trans_seqs = {}
    for chr in gene_ann.__transcripts_ann:
        seq = []
        for transcript in gene_ann.__transcripts_ann[chr]:  
            seq.append(refs.fetch(transcript.chr, transcript.start, transcript.end))        
                
        trans_seqs.setdefault(transcript, "".join(seq))
    return trans_seqs

#split gene order string to get gene names
# MYH8_cds>>MYH13_intron;RP11-401O9.3_utr5>>AK097500_utr3;CTC-297N7.11_utr5;RP11-799N11.1_utr5
def get_gene_names_from_gene_order(gene_order_string):
    # gene name list
    g_driver = []
    g_passenger = []
    tmp = gene_order_string.split(";")
    for t in tmp:
        tmp2 = t.split(">>")
        for n in range(0, len(tmp2)):
            t2 = tmp2[n]
            if n == 0:
                g = g_driver
            else:
                g = g_passenger
            tmp3 = t2.split(",")
            for t3 in tmp3:
                if "_" in t3:
                    tmp4 = t3.split("_")
                    g.append(tmp4[0])
    return g_driver, g_passenger

# build fusion sequence, with up/downstream genome sequence, splicing not considered by far     
def get_fusion_seq(fusion, ref, seg_len):

    chr1 = fusion.chr1
    chr2 = fusion.chr2
    bp1 = fusion.pos1
    bp2 = fusion.pos2
    strand1 = fusion.strand1
    strand2 = fusion.strand2

    refs = pysam.FastaFile(ref)
    if not chr1 in refs:
        if chr1.startswith("chr"):
            chr1 = chr1[3:]
            chr2 = chr2[3:]
        else:
            chr1 = "chr" + chr1
            chr2 = "chr" + chr2
    if not chr1 in refs:
        print >> sys.stderr, "unmatched reference."
        sys.exit(1)


    if strand1 == "+": ## ->|
        win_start1 = bp1 -seg_len
        win_end1 = bp1
        seq1 = refs.fetch(chr1, win_start1, win_end1)
    else: ## |<-
        win_start1 = bp1 - 1
        win_end1 = bp1 +seg_len - 1
        seq1 = sequtils.rc_seq(refs.fetch(chr1, win_start1, win_end1), "rc")

    if strand2 == "+": ## <-|
        win_start2 = bp2 - seg_len
        win_end2 = bp2
        seq2 = sequtils.rc_seq(refs.fetch(chr2, win_start2, win_end2), "rc")
    else: ## |->
        win_start2 = bp2 - 1
        win_end2 = bp2 + seg_len - 1
        seq2 = refs.fetch(chr2, win_start2, win_end2)
    fusion.seq1 = seq1.upper()
    fusion.seq2 = seq2.upper()
    if not fusion.seq1:
        fusion.seq1 = "NA"
    if not fusion.seq2:
        fusion.seq2 = "NA"
        
    refs.close()
    #return seq1.upper(), seq2.upper()

def format_chr(chr0, ref):  
    formatted_chr = chr0
    if not formatted_chr in ref:
        if formatted_chr.startswith("chr"):
            formatted_chr = formatted_chr[3:]
        else:
            formatted_chr = "chr" + formatted_chr
    
    return formatted_chr

def fetch_noncoding_seqs(gene_ann, cff_fusion, ref, bp, gene_order):
    """
    Extracts sequences surounding breakpoint for genes which are not in annotation ("NA", or otherwise)
    Both forward and reverse sequences are extracted, and both will be used to query transcriptome,
    since strand information is not available for this gene 
    """
    seqs = []
    if gene_order == "head":
        # extract seqs from lhs and rhs of breakpoint 
        left_seq = ref.fetch(format_chr(cff_fusion.chr1, ref), bp-101, bp-1) 
        right_seq = ref.fetch(format_chr(cff_fusion.chr1, ref), bp-1, bp+99)
        #generate fwd seqs
        fwd_fusion_seq = left_seq
        trans_seq1, trans_seq2 = left_seq, right_seq #+"head_fwd" 
        seqs.append((trans_seq1, trans_seq2, fwd_fusion_seq))
        #generate rev seqs
        rev_fusion_seq = sequtils.rc_seq(right_seq, "rc")
        trans_seq1, trans_seq2 = sequtils.rc_seq(right_seq, "rc"), sequtils.rc_seq(left_seq, "rc") #+ "head_rev"
        seqs.append((trans_seq1, trans_seq2, rev_fusion_seq))

    elif gene_order == "tail":
        left_seq = ref.fetch(format_chr(cff_fusion.chr2, ref), bp-101, bp-1)
        right_seq = ref.fetch(format_chr(cff_fusion.chr2, ref), bp-1, bp+99)
        #generate fwd seqs
        fwd_fusion_seq = right_seq # + "tail_fwd" 
        trans_seq1, trans_seq2 = left_seq, right_seq
        #trans_seq1, trans_seq2 =  
        seqs.append((trans_seq1, trans_seq2, fwd_fusion_seq))
        #generate rev seqs
        rev_fusion_seq = sequtils.rc_seq(left_seq, "rc") #+"tail_rev"
        trans_seq1, trans_seq2 = sequtils.rc_seq(right_seq, "rc"), sequtils.rc_seq(left_seq, "rc")
        #print(trans_seq1, trans_seq2)
        #exit(0)
        seqs.append((trans_seq1, trans_seq2, rev_fusion_seq))
    return seqs 

# for a given gene, according to gene_order and bp get all its transcript sequences and 100bp potential fusion sequences (only for head gene)
def build_transcript_and_fusion_seq(gene_ann, cff_fusion, ref, bp, gene_order):
    # get gene name from cff_fusion object
    if gene_order == "head":
        gene_name = cff_fusion.reann_gene1 
    elif gene_order == "tail":
        gene_name = cff_fusion.reann_gene2

    #noncoding gene
    noncoding=0
    # get all transcripts for head gene 
    seqs = []
    gene_intervals = gene_ann.get_gene_interval(gene_name)
    tids = gene_intervals.transcript_ids
    #print(tids)
    if not tids:
        print >> sys.stderr, """Gene not in annotation: {}. """.format(gene_name) + "Building reference sequence from noncoding region"
        noncoding=1
        seqs = fetch_noncoding_seqs(gene_ann, cff_fusion, ref, bp, gene_order)
        return seqs
        
    is_head = (gene_order == "head")
    for transcript_id in tids:
        trans_ann_list = gene_ann.get_transcripts_ann(transcript_id) # all GeneBed annotations of transcript_id
        if not trans_ann_list[0].start <= bp <= trans_ann_list[-1].end: #skip the transcript that doesnt include the breakpoint
            continue
        is_fw = (trans_ann_list[0].strand == "f")
        trans_seqs = []
        fusion_seqs = []
        if is_head:
            if is_fw: # head gene on forward strand
                for ann in trans_ann_list:
                    # breakpoint in current bed annotation
                    if ann.start <= bp <= ann.end:
                        seq = ref.fetch(format_chr(ann.chr, ref), ann.start-1, bp)
                        fusion_seqs.append(seq)
                        if ann.type == "intron": # assume a intron-retention
                            trans_seqs.append(seq)
                    if ann.type != "intron":
                        seq = ref.fetch(format_chr(ann.chr, ref), ann.start-1, ann.end)
                        trans_seqs.append(seq)
                        if ann.end < bp:
                            fusion_seqs.append(seq)
            else: # head gene on reverse strand
                for ann in trans_ann_list:
                    # breakpoint in current bed annotation
                    if ann.start <= bp <= ann.end:
                        seq = ref.fetch(format_chr(ann.chr, ref), bp-1, ann.end)
                        fusion_seqs.append(seq)
                        if ann.type == "intron": # breakpoint in intron, cant build seq
                            trans_seqs.append(seq)
                        
                    if ann.type != "intron": # for now do not consider intonic breakpoint
                        seq = ref.fetch(format_chr(ann.chr, ref), ann.start-1, ann.end)
                        trans_seqs.append(seq)
                        if ann.start > bp:
                            fusion_seqs.append(seq)

        else: # tail gene, do not build transcript seq
            if is_fw: # tail gene on forward strand
                for ann in trans_ann_list:
                    # breakpoint in current bed annotation
                    if ann.start <= bp <= ann.end:
                        seq = ref.fetch(format_chr(ann.chr, ref), bp-1, ann.end)
                        fusion_seqs.append(seq)
                    elif ann.start > bp and ann.type != "intron": # for now do not consider intonic breakpoint
                        seq = ref.fetch(format_chr(ann.chr, ref), ann.start-1, ann.end)
                        fusion_seqs.append(seq)
                        
            else:
                for ann in trans_ann_list:
                    # breakpoint in current bed annotation
                    if ann.start <= bp <= ann.end:
                        seq = ref.fetch(format_chr(ann.chr, ref), ann.start-1, bp)
                        fusion_seqs.append(seq)
                        break
                    elif ann.type != "intron": # for now do not consider intonic breakpoint
                        seq = ref.fetch(format_chr(ann.chr, ref), ann.start-1, ann.end)
                        fusion_seqs.append(seq)

        trans_seq = ("".join(trans_seqs)).upper()
        fusion_seq = ("".join(fusion_seqs)).upper()
            
        if fusion_seq:
            seq_len = 100
            if not is_fw:
                trans_seq = sequtils.rc_seq(trans_seq, "rc")
                fusion_seq = sequtils.rc_seq(fusion_seq, "rc")

            if is_head:
                fusion_bp = len(fusion_seq)
                fusion_seq = fusion_seq[-seq_len:]
                trans_seq1 = trans_seq[max(0, fusion_bp-seq_len):fusion_bp]
                trans_seq2 = trans_seq[fusion_bp:fusion_bp+seq_len]

            else:
                fusion_bp = len(trans_seq) - len(fusion_seq)
                fusion_seq = fusion_seq[:seq_len]
                trans_seq1 = trans_seq[max(0, fusion_bp-seq_len):fusion_bp]
                trans_seq2 = trans_seq[fusion_bp:fusion_bp+seq_len]

            seqs.append((trans_seq1, trans_seq2, fusion_seq))
    return seqs

# for given gene annotation bed file, build junctions sequences for echo transcript
# ref is pysam.FastaFile

#chr1    24683494        24685032        ENST00000003583 utr3    0       r       STPG1   ENSG00000001460
#chr1    24685032        24685109        ENST00000003583 cds     0       r       STPG1   ENSG00000001460
#chr1    24685109        24687340        ENST00000003583 intron  0       r       STPG1   ENSG00000001460
def output_trans_seq(trans_seqs, bp_list, cur_trans_id, cur_trans_strand, junc_seq_len):
    trans_seq = "".join(trans_seqs)
    junc_id = 0
    for bp in bp_list:
        junc_id += 1
        seq1 = trans_seq[max(0, bp-junc_seq_len):bp]
        seq2 = trans_seq[bp:min(len(trans_seq), bp+junc_seq_len)]
        seq = seq1 + seq2
        if cur_trans_strand == "r":
            seq = sequtils.rc_seq(seq, "rc")
        ref_name = ">" + cur_trans_id + "_" + str(junc_id)
        print ref_name
        print seq

    
def build_junction_seq_for_gene_bed(ref, gene_bed):
    cur_trans_id = ""
    cur_trans_strand = ""
    cur_trans_chr = ""
    junc_seq_len = 100
    trans_seqs = []
    bp_list = []
    for line in open(gene_bed, "r"):
        ann = GeneBed(line)
        if not format_chr(ann.chr, ref) in ref:
            continue
        if cur_trans_id != ann.transcript_id:
            #output transcript junction sequences
            output_trans_seq(trans_seqs, bp_list, cur_trans_id, cur_trans_strand, junc_seq_len)

            cur_trans_id = ann.transcript_id
            cur_trans_strand = ann.strand
            cur_trans_chr = ann.chr
            bp_list = []
            trans_seqs = []
        if ann.type != "intron":
            seq = ref.fetch(format_chr(ann.chr, ref), ann.start-1, ann.end) 
            trans_seqs.append(seq)
        else:
            bp_list.append(len("".join(trans_seqs)))
    output_trans_seq(trans_seqs, bp_list, cur_trans_id, cur_trans_strand, junc_seq_len)
            
                