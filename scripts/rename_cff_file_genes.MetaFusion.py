#!/usr/bin/env python
import os
import sys
import pygeneann_MetaFusion as pygeneann
import pandas as pd
import sequtils
import argparse
import re

parser = argparse.ArgumentParser()
parser.add_argument('cff_file', action='store', help='CFF file before annotation. if there are multiple gene names in a field, names MUST be comma-separated lists')
parser.add_argument('gene_info_file', action='store', help='Homo_sapiens.gene_info')

args = parser.parse_args()
cff_file = args.cff_file
gene_info_file=args.gene_info_file

# tab-delimited file with a column Symbol == HGNC_Symbol and Synonyms== ENSG_id

#Open Gene info file and create df
ncbi_gene_info_file = open(gene_info_file)
df = pd.read_csv(ncbi_gene_info_file, sep='\t')

#Create one row for each alias
df = pd.concat([pd.Series(row['Symbol'], row['Synonyms'].split('|'))
	for _, row in df.iterrows()]).reset_index()
df.rename(columns = {list(df)[0]: 'Alias', list(df)[1]: 'HGNC_Symbol'}, inplace = True)

# Remove bracketed gene information from arriba
def clean_weird_input(gene_name): 
    #$ cat $cluster | sed 's/(.\+)//g' | sed 's/\//,/g' 
    # remove arriba brackets, e.g. RNU6-121P(14818),AQP10(8524)
    gene_name = re.sub('\([0-9]+\)', '', gene_name) 
    # remove integrate slash-delimited paralogs
    gene_name = re.sub('/', ',', gene_name) 
    return gene_name 

def alias2hgnc(df, query):
    #check to see if query is already HGNC symbol
    hgnc = df.loc[df.HGNC_Symbol == query].HGNC_Symbol.values.tolist()
    is_hgnc = True if len(hgnc) > 0 else False 
    if is_hgnc: 
        return query
    else:
        hgnc_lst = df.loc[df.Alias == query].HGNC_Symbol.values.tolist()
        return "NA" if len(hgnc_lst) == 0 else hgnc_lst[0]

# In arriba, can have more than 1 gene on fusion. This will select the "more likely" Symbol
def select_gene_name(gene_lst):
    """takes a list of gene names and selects a single name"""
    #handle list of genes from caller
    gene_hgnc_lst = []
    for gene in gene_lst:
        hgnc_gene = alias2hgnc(df, gene)
        gene_hgnc_lst.append(hgnc_gene)
    try: return [gen for gen in gene_hgnc_lst if gen != "NA"][0]
    except: return ','.join([gen for gen in gene_hgnc_lst if gen != "NA"])

for line in open(cff_file, "r"):
    fusion = pygeneann.CffFusion(line)
    #clean weird input
    fusion.t_gene1 = clean_weird_input(fusion.t_gene1) 
    fusion.t_gene2 = clean_weird_input(fusion.t_gene2) 
    # t_gene1
    t_gene1 = select_gene_name(fusion.t_gene1.split(","))
    ### If no gene name is selected, return the callers gene name
    # This should not occur if you use the generate gene info file with all gtfs from input callers
    if len(t_gene1) != 0: fusion.t_gene1 = t_gene1 
    # t_gene2
    t_gene2 = select_gene_name(fusion.t_gene2.split(","))
    if len(t_gene2) != 0: fusion.t_gene2 = t_gene2 

    print(fusion.tostring())    
    #gene2_lst = fusion.t_gene2.split(",")
    #hgnc_gene1 = alias2hgnc(df, fusion.t_gene1) 
    #hgnc_gene2 = alias2hgnc(df, fusion.t_gene2) 
    #print( hgnc_gene1,hgnc_gene2)
    