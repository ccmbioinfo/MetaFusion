#!/usr/bin/env python
import sys
import pandas as pd
import numpy as np
sys.path.append("/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin")
import  pygeneann_reads_capture_DEV as pygeneann
import pybedtools.bedtool as bedtools
import itertools
import sequtils
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('cff', action='store', help='CFF file')

args = parser.parse_args()
cff=args.cff
#cff="/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/testing_pipeline/NEG_CONTROL_BEERS.benchmark.June-5-2020/merged.cff.renamed.reann"

#INTERSECT FUSIONS BY BREAKPOINTS
def intersect_fusions_by_breakpoints():
    lines=[line for line in open(cff, "r")]
    fusion=pygeneann.CffFusion(lines[0])
    header=fusion.zone1_attrs + fusion.zone2_attrs + fusion.zone3_attrs + fusion.zone4_attrs
    df_cff=pd.read_csv(cff, sep='\t', keep_default_na=False, index_col=False, names=header)
    
    
    #create BedTools object with appropriate column names
    print >> sys.stderr, "create BedTools object with appropriate column names"
    df_bed=df_cff[['chr1','pos1','pos1','chr2','pos2','pos2', 'fusion_id']]
    df_bed.columns=['chr1','pos1','pos1_2','chr2','pos2','pos2_2', 'fusion_id']
    df_bed.loc[:,['pos1_2','pos2_2']] +=1
    df_bed=bedtools.BedTool.from_dataframe(df_bed)
    
    #Intersect fusions: NOTE: only keeps fusions that intersect
    #print >> sys.stderr, "Intersect fusions: NOTE: only keeps fusions that intersect"
    #df_intersect=df_bed.pair_to_pair(df_bed, slop=100, rdn=True)
    print >> sys.stderr, "Intersect fusions: NOTE: rdn=False, keeps self-intersections"
    df_intersect=df_bed.pair_to_pair(df_bed, slop=100, rdn=False)
    df=df_intersect.to_dataframe(header=None).iloc[:,0:14]
    df.columns = ['chr1','pos1','pos1_2','chr2','pos2','pos2_2', 'fusion_id', 'chr1_1','pos1_1','pos1_2_1','chr2_1','pos2_1','pos2_2_1', 'fusion_id_lst'] 
    df=df[['fusion_id','fusion_id_lst']]
    #write paired F_IDs to tsv
    return df

df = intersect_fusions_by_breakpoints()
df.to_csv(sys.stdout,header=True,index=True, sep="\t")
#exit(0)
#df.to_csv("/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/test_R_graph_cluster/" + sample_name + ".FIDs.intersections.tsv", sep = "\t")

# CLUSTER GENES

def intersect_fusions_by_genes(cff_file):
    fusion_dict = {}
    fusion_list_for_bp_cmp = []
    common_key_dict = {}
    # cluster fusions by gene pairs, save in fusion_dict
    for line in open(cff_file, "r"):
        if line.startswith("#"):
            continue
        fusion = pygeneann.CffFusion(line)
        if fusion.t_gene1 == "NA" or fusion.t_gene2 == "NA":
            continue
        else:
            key = ",".join(sorted([fusion.t_gene1 + "|" + fusion.chr1, fusion.t_gene2+ "|" + fusion.chr2])) 
            fusion_dict.setdefault(key, []).append(fusion.fusion_id)
    return fusion_dict

fusion_dict = intersect_fusions_by_genes(cff)
#print([(key, len(fusion_dict[key])) for key in fusion_dict.keys()])

#exit(0)
count = df.shape[0] + 1 
for key in fusion_dict.keys():
    lst=fusion_dict[key]
    edges=list(itertools.permutations(lst, 2))
    for edge in edges:
        print("\t".join([str(count)] + list(edge)))
        count += 1

# RUN GRAPH CLUSTERING
#/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/test_R_graph_cluster/DIPG.all_edges.tsv

# OUTPUT CLUSTERED FUSIONS


#self.output_clustered_fusions(fusion_list, "Gene_Cluster")
#print(count)
