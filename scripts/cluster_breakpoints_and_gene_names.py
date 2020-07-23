import pandas as pd
import numpy as np
import  pygeneann_reads_capture_DEV as pygeneann
import pybedtools.bedtool as bedtools
import itertools


#cff="/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/testing_pipeline/outdir.TEST_FUSIONS.May-19-2020-WITH_RENAME/BCL3--CTB-171A8.1.cff.renamed.reann"
cff="/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/testing_pipeline/outdir.TEST_FUSIONS.May-19-2020-WITH_RENAME/BCL3--CTB-171A8.1.cff.renamed.reann.MOD"
#cff="/hpf/largeprojects/ccmbio/mapostolides/DIPG/run_DIPG_samples/DIPG_output_SEPT4-2019_TESTING_readset_fastq_start/fusions/cff/merged.cff.reann.dnasupp"
cff="/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/testing_pipeline/DIPG.May-26-2020_genebed_renamed.ENSG_gene_bed/merged.cff.renamed.reann"
sample_name="DIPG"
#cff="/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/testing_pipeline/SIM45.SIM52.benchmark.May-25-2020_genebed_renamed.ENSG/merged.cff.renamed.reann"
#sample_name="SIM50"

lines=[line for line in open(cff, "r")]
#fusion=pygeneann.CffFusion(lines[0])
#header=fusion.zone1_attrs + fusion.zone2_attrs + fusion.zone3_attrs + fusion.zone4_attrs
#df_cff=pd.read_csv(cff, sep='\t', keep_default_na=False, index_col=False, names=header)
#
##INTERSECT FUSIONS BY BREAKPOINTS
#
##create BedTools object with appropriate column names
#print("create BedTools object with appropriate column names")
#df_bed=df_cff[['chr1','pos1','pos1','chr2','pos2','pos2', 'fusion_id']]
#df_bed.columns=['chr1','pos1','pos1_2','chr2','pos2','pos2_2', 'fusion_id']
#df_bed.loc[:,['pos1_2','pos2_2']] +=1
#df_bed=bedtools.BedTool.from_dataframe(df_bed)
#
##Intersect fusions: NOTE: only keeps fusions that intersect
#print("Intersect fusions: NOTE: only keeps fusions that intersect")
#df_intersect=df_bed.pair_to_pair(df_bed, slop=100, rdn=True)
#df=df_intersect.to_dataframe(header=None).iloc[:,0:14]
#df.columns = ['chr1','pos1','pos1_2','chr2','pos2','pos2_2', 'fusion_id', 'chr1_1','pos1_1','pos1_2_1','chr2_1','pos2_1','pos2_2_1', 'fusion_id_lst'] 
#df=df[['fusion_id','fusion_id_lst']]
##write paired F_IDs to tsv
#df.to_csv("/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/test_R_graph_cluster/" + sample_name + ".FIDs.intersections.tsv", sep = "\t")


# CLUSTER GENES

def cluster_fusions_by_genes(cff_file):
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
    # output clustered fusions
    #for key in fusion_dict:
    #    fusion_list = fusion_dict[key]

        #self.output_clustered_fusions(fusion_list, "Gene_Cluster")
    return fusion_dict

fusion_dict = cluster_fusions_by_genes(cff)
#print([(key, len(fusion_dict[key])) for key in fusion_dict.keys()])
count = 351896 
for key in fusion_dict.keys():
    #lst=fusion_dict['CELSR1|chr22,SPECC1L|chr22']
    lst=fusion_dict[key]
    #print(list(itertools.permutations(lst, 2)))
    #count += len(list(itertools.permutations(lst, 2)))
    edges=list(itertools.permutations(lst, 2))
    for edge in edges:
        print("\t".join([str(count)] + list(edge)))
        count += 1
    
#print(count)
