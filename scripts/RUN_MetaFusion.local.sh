#!/bin/bash
#NOTES
#source /hpf/largeprojects/ccmbio/mapostolides/MODULES/miniconda3/etc/profile.d/conda.sh
#conda activate metafusion

#REFERENCE FILES FILES
test_dir=/Users/mapostolides/MetaFusion/RUNS
gene_bed=/Users/mapostolides/MetaFusion/reference_files/ens_known_genes.renamed.ENSG.bed
gene_info=/Users/mapostolides/MetaFusion/reference_files/Homo_sapiens.gene_info
genome_fasta=/Users/mapostolides/MetaFusion/reference_files/human_g1k_v37_decoy.fasta
recurrent_bedpe=/Users/mapostolides/MetaFusion/reference_files/blacklist_breakpoints.bedpe


date=July-23-2020
#DATASETS
sim45_sim52=1

# SIM45.SIM52.combined
if [ $sim45_sim52 -eq 1 ]; then
echo SIM45.SIM52
outdir=$test_dir/SIM45.SIM52.benchmark.$date.MetaFusion
cff=/Users/mapostolides/MetaFusion/RUNS/cff_files/sim45.sim52.merged.cff
truth_fusions=sim45.sim52.truth_set.dat

sh MetaFusion.sh --outdir $outdir --cff $cff  --gene_bed $gene_bed  --gene_info $gene_info --truth_set $truth_fusions  --num_tools=2  --genome_fasta $genome_fasta --recurrent_bedpe $recurrent_bedpe
num_tools=2
#./MetaFusion.sh $outdir $cff $gene_bed $truth_fusions $num_tools $gene_info $genome_fasta
#sh MetaFusion.sh --outdir $outdir \
#                --cff $cff \
#                --gene_bed $gene_bed \
#                --gene_info $gene_info \
#                --truth_set $truth_fusions \
#                --num_tools=2 \
#                --genome_fasta $genome_fasta
fi

