#!/bin/bash

#Change date to current date. Can also add tag to this string for multiple runs
date=Sept-21-2020-Singularity-V3
topdir=/hpf/largeprojects/ccmbio/mapostolides
#DATASETS
sim45_sim52=1
#brca_4=1
#beers_neg=1
#sim_50=1
#sim101=1
#melanoma=1

fusiontools=$topdir/MetaFusion/scripts
#REFERENCE FILES FILES
runs_dir=$topdir/MetaFusion/RUNS
mkdir $runs_dir
gene_bed=$topdir/MetaFusion/reference_files/ens_known_genes.renamed.ENSG.bed
gene_info=$topdir/MetaFusion/reference_files/Homo_sapiens.gene_info
genome_fasta=$topdir/MetaFusion/reference_files/human_g1k_v37_decoy.fasta
recurrent_bedpe=$topdir/MetaFusion/reference_files/blocklist_breakpoints.bedpe


# SIM45.SIM52.combined
if [ $sim45_sim52 -eq 1 ]; then
echo SIM45.SIM52
outdir=$runs_dir/SIM45.SIM52.benchmark.$date
echo generating output in $outdir
mkdir $outdir
cff=$topdir/MetaFusion/test_data/cff/dream.sim45.sim52.cff
truth_fusions=$topdir/MetaFusion/test_data/truth_sets/dream.sim45.sim52.truth_set.dat

bash MetaFusion.sh --outdir $outdir \
                 --cff $cff  \
                 --gene_bed $gene_bed \
                 --fusion_annotator \
                 --genome_fasta $genome_fasta \
                 --gene_info $gene_info \
                 --truth_set $truth_fusions \
                 --num_tools=2  \
                 --recurrent_bedpe $recurrent_bedpe \
                 --scripts $fusiontools 
fi

