#!/usr/bin/env python
import sys
import pandas as pd
import numpy as np
import pygeneann_MetaFusion as pygeneann
import pybedtools.bedtool as bedtools
import itertools
import sequtils
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('cff', action='store', help='CFF file')

args = parser.parse_args()
cff=args.cff

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
    print >> sys.stderr, "Intersect fusions: NOTE: rdn=False, keeps self-intersections"
    df_intersect=df_bed.pair_to_pair(df_bed, slop=100, rdn=False)
    df=df_intersect.to_dataframe(header=None).iloc[:,0:14]
    df.columns = ['chr1','pos1','pos1_2','chr2','pos2','pos2_2', 'fusion_id', 'chr1_1','pos1_1','pos1_2_1','chr2_1','pos2_1','pos2_2_1', 'fusion_id_lst'] 
    df=df[['fusion_id','fusion_id_lst']]
    #write paired F_IDs to tsv
    return df

#df = intersect_fusions_by_breakpoints()
#df.to_csv(sys.stdout,header=True,index=True, sep="\t")

# CLUSTER EXONS 

def intersect_fusions_by_exons(cff):
    fusion_dict = {}
    # cluster fusions by gene pairs, save in fusion_dict
    for line in open(cff, "r"):
        if line.startswith("#"):
            continue
        fusion = pygeneann.CffFusion(line)
        # skip those fusions which lack "closest exons" near at least one of their breakpoints
        if fusion.closest_exon1 == "NA" or fusion.closest_exon2 == "NA":
            continue
        else:
            key = ",".join(sorted([fusion.closest_exon1, fusion.closest_exon2])) 
            # append fusion_id twice so it intersects with itself
            fusion_dict.setdefault(key, []).append(fusion.fusion_id)
            fusion_dict.setdefault(key, []).append(fusion.fusion_id)
    return fusion_dict

fusion_dict = intersect_fusions_by_exons(cff)
#print(fusion_dict)
#exit(0)
#count = df.shape[0] + 1 
count = 0
#print header
print("\t" + "fusion_id" + "\t" + "fusion_id_lst")
for key in fusion_dict.keys():
    lst=fusion_dict[key]
    #print(lst)
    edges=list(itertools.permutations(lst, 2))
    #print(edges)
    for edge in edges:
        print("\t".join([str(count)] + list(edge)))
        count += 1
