#!/bin/bash

#Change date to current date. Can also add tag to this string for multiple runs
date=$1
topdir=$2

fusiontools=$topdir/MetaFusion/scripts
#REFERENCE FILES FILES
runs_dir=$topdir/MetaFusion/RUNS
gene_bed=$topdir/MetaFusion/reference_files/ens_known_genes.renamed.ENSG.bed
gene_info=$topdir/MetaFusion/reference_files/Homo_sapiens.gene_info
genome_fasta=$topdir/MetaFusion/reference_files/human_g1k_v37_decoy.fasta
recurrent_bedpe=$topdir/MetaFusion/reference_files/blocklist_breakpoints.bedpe

outdir=$3
echo generating output in $outdir
mkdir $outdir
cff=$4
truth_fusions=$5

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

