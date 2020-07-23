#!/bin/bash
module load R/3.5.1
#module load python/2.7.11
fusiontools_dir=/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin
cff=$1
outdir=$2

#Generate intersections file for both gene name and breakpoint
fid_intersection_file=$outdir/FID.intersections.tsv
python $fusiontools_dir/intersect_breakpoints_and_gene_names.py $cff > $fid_intersection_file 

# Graph clustering
fid_clusters_file=$outdir/FID.clusters.tsv
Rscript $fusiontools_dir/cluster_intersections.R $fid_intersection_file $fid_clusters_file 

# generate cluster file using clustered FIDs and cff file
python $fusiontools_dir/generate_cluster_file.py $cff $fid_clusters_file
