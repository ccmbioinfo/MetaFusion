#!/usr/bin/env python
#/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/validation_pipeline/Pipeline-scripts/rename_cff_file_genes-GENAP.py
import os
import sys
import pygeneann_reads_capture_DEV as pygeneann
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

# tab-delimited file with a column contating all aliases for each gene
#tax_id	GeneID	Symbol	LocusTag	Synonyms	dbXrefs	chromosome	map_location	description	type_of_gene	Symbol_from_nomenclature_authority	Full_name_from_nomenclature_authority	Nomenclature_status	Other_designations	Modification_date	Feature_type
#9606	1	A1BG	-	A1B|ABG|GAB|HYST2477	MIM:138670|HGNC:HGNC:5|Ensembl:ENSG00000121410	19	19q13.43	alpha-1-B glycoprotein	protein-coding	A1BG	alpha-1-B glycoprotein	O	alpha-1B-glycoprotein|HEL-S-163pA|epididymis secretory sperm binding protein Li 163pA	20200313	-

#Open NCBI file and create df
#ncbi_gene_info_file = open("/hpf/largeprojects/ccmbio/mapostolides/gene_fusion/pipeline/config_reference_files/Homo_sapiens.gene_info") 
ncbi_gene_info_file = open(gene_info_file)
df = pd.read_csv(ncbi_gene_info_file, sep='\t')

#https://stackoverflow.com/questions/26336251/pandas-rename-single-dataframe-column-without-knowing-column-name
#Create one row for each alias
df = pd.concat([pd.Series(row['Symbol'], row['Synonyms'].split('|'))
	for _, row in df.iterrows()]).reset_index()
df.rename(columns = {list(df)[0]: 'Alias', list(df)[1]: 'HGNC_Symbol'}, inplace = True)

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

def select_gene_name(gene_lst):
    """takes a list of gene names and selects a single name"""
    #handle list of genes from caller
    gene_hgnc_lst = []
    for gene in gene_lst:
        hgnc_gene = alias2hgnc(df, gene)
        gene_hgnc_lst.append(hgnc_gene)
    try: return [gen for gen in gene_hgnc_lst if gen != "NA"][0]
    except: return ','.join([gen for gen in gene_hgnc_lst if gen != "NA"])

#cff_file ="/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/reann_cff_fusion_testing/FOSB--AADACL2/FOSB--AADACL2.cff.MOD"
for line in open(cff_file, "r"):
    fusion = pygeneann.CffFusion(line)
    #clean weird input
    fusion.t_gene1 = clean_weird_input(fusion.t_gene1) 
    fusion.t_gene2 = clean_weird_input(fusion.t_gene2) 
    # t_gene1
    t_gene1 = select_gene_name(fusion.t_gene1.split(","))
    if len(t_gene1) != 0: fusion.t_gene1 = t_gene1 
    # t_gene2
    t_gene2 = select_gene_name(fusion.t_gene2.split(","))
    if len(t_gene2) != 0: fusion.t_gene2 = t_gene2 

    print(fusion.tostring())    
    #gene2_lst = fusion.t_gene2.split(",")
    #hgnc_gene1 = alias2hgnc(df, fusion.t_gene1) 
    #hgnc_gene2 = alias2hgnc(df, fusion.t_gene2) 
    #print( hgnc_gene1,hgnc_gene2)
