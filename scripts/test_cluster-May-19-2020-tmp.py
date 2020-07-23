import pandas as pd
import numpy as np
import  pygeneann_reads_capture_DEV as pygeneann
import pybedtools.bedtool as bedtools
#cff="/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/testing_pipeline/outdir.TEST_FUSIONS.May-19-2020-WITH_RENAME/BCL3--CTB-171A8.1.cff.renamed.reann"
sample_name="TEST"
cff="/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/testing_pipeline/outdir.TEST_FUSIONS.May-19-2020-WITH_RENAME/BCL3--CTB-171A8.1.cff.renamed.reann.MOD"
#cff="/hpf/largeprojects/ccmbio/mapostolides/DIPG/run_DIPG_samples/DIPG_output_SEPT4-2019_TESTING_readset_fastq_start/fusions/cff/merged.cff.reann.dnasupp"
#sample_name="DIPG"
#cff="/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/testing_pipeline/SIM45.SIM52.benchmark.May-25-2020_genebed_renamed.ENSG/merged.cff.renamed.reann"
#sample_name="SIM50"

lines=[line for line in open(cff, "r")]
fusion=pygeneann.CffFusion(lines[0])
header=fusion.zone1_attrs + fusion.zone2_attrs + fusion.zone3_attrs + fusion.zone4_attrs
df_cff=pd.read_csv(cff, sep='\t', keep_default_na=False, index_col=False, names=header)

#INTERSECT FUSIONS BY BREAKPOINTS

#create BedTools object with appropriate column names
print("create BedTools object with appropriate column names")
df_bed=df_cff[['chr1','pos1','pos1','chr2','pos2','pos2', 'fusion_id']]
df_bed.columns=['chr1','pos1','pos1_2','chr2','pos2','pos2_2', 'fusion_id']
df_bed.loc[:,['pos1_2','pos2_2']] +=1
df_bed=bedtools.BedTool.from_dataframe(df_bed)

#Intersect fusions: NOTE: only keeps fusions that intersect
print("Intersect fusions: NOTE: only keeps fusions that intersect")
df_intersect=df_bed.pair_to_pair(df_bed, slop=100, rdn=True)
df=df_intersect.to_dataframe(header=None).iloc[:,0:14]
df.columns = ['chr1','pos1','pos1_2','chr2','pos2','pos2_2', 'fusion_id', 'chr1_1','pos1_1','pos1_2_1','chr2_1','pos2_1','pos2_2_1', 'fusion_id_lst'] 
df=df[['fusion_id','fusion_id_lst']]
#write paired F_IDs to tsv
df.to_csv("/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/test_R_graph_cluster/" + sample_name + ".FIDs.intersections.tsv", sep = "\t")

#print(df)
exit(0)
#remove duplicate rows before grouping (i.e. A-->B and B-->A are duplicates)
print("remove duplicate rows before grouping (i.e. A-->B and B-->A are duplicates)")
df=pd.DataFrame(np.sort(df.values, axis=1), columns=df.columns).drop_duplicates()
#print(df)

new = pd.DataFrame(columns=df.columns)
#ids_so_far=[]
for i in df.index: 
  try: df = df[df.fusion_id != df.loc[i].fusion_id_lst]
  except: continue
  #if df.loc[i].fusion_id in ids_so_far: continue
  #ids_so_far.append(df.loc[i].fusion_id_lst)
#  print(df.loc[i])
  #new=new.append(df.loc[i])
#df=new
##group all intersecting fusions
print("group all intersecting fusions")
grouped=df.groupby('fusion_id')['fusion_id_lst'].apply(lambda x: ','.join(x))
df=pd.DataFrame(grouped)
#print(df)

#convert "fusion_id" rownames to column
df.reset_index(inplace=True)
##join query f_id with intersecting list
df["fusion_id_lst"] = df["fusion_id"] + "," + df["fusion_id_lst"] 

print(df)
df.to_csv("/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/test_R_graph_cluster/DIPG.FIDs.clustered.tsv", sep = "\t")
print("done")
#print("getting unique fusion_id lists")
#fid_list=df.fusion_id_lst.to_list()
#print(len(fid_list))
#fid_list=[fid for fid in fid_list if len(filter(lambda x: fid in x,fid_list)) == 1]
#print(len(fid_list))

#df=pd.DataFrame({'fusion_id_lst': fid_list})
#print(df)
##exit(0)
## get fusion_ids matching gene names
#print(" get fusion_ids matching gene names")
#for i in range(0, len(df.index)):
#  fusion_id_lst=df.iloc[i][0].split(",")
#  df_tmp=df_cff[df_cff.fusion_id.isin(fusion_id_lst)]
#  head_genes=list(set(df_tmp['t_gene1'].tolist()))
#  tail_genes=list(set(df_tmp['t_gene2'].tolist()))
#  #check gene name matches for both orientations
#  for j in range(0, len(df_cff.index)): 
#    fusion=df_cff.iloc[j]
#    if fusion.t_gene1 in head_genes and fusion.t_gene2 in tail_genes: 
#      fusion_id_lst.append(fusion.fusion_id)
#    elif fusion.t_gene2 in head_genes and fusion.t_gene1 in tail_genes:
#      fusion_id_lst.append(fusion.fusion_id)
#  print(set(fusion_id_lst))
#compare set of names for HEAD and TAIL genes 

#SCRAP
#execfile('test_cluster-May-19-2020-tmp.py')
#print(df_grouped.iloc[:,6])
#print(df.to_string())

#header=['chr1', 'pos1', 'strand1', 'chr2', 'pos2', 'strand2', 'library', 'sample_name', 'sample_type', 'disease', 'tool', 'split_cnt', 'span_cnt', 't_gene1', 't_area1', 't_gene2', 't_area2', 'category', 'reann_gene1', 'reann_type1', 'reann_gene2', 'reann_type2', 'gene1_on_bdry', 'gene1_close_to_bndry', 'gene2_on_bdry', 'gene2_close_to_bndry', 'score', 'coding_id_distance', 'gene_interval_distance', 'dna_support', 'fusion_id', 'seq1', 'seq2', 'is_inframe', 'splice_site1', 'splice_site2', 'captured_reads']
