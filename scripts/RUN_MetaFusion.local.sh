#!/bin/bash
#NOTES
#source /hpf/largeprojects/ccmbio/mapostolides/MODULES/miniconda3/etc/profile.d/conda.sh
#conda activate metafusion

fusiontools=/Users/mapostolides/MetaFusion/scripts
#REFERENCE FILES FILES
test_dir=/Users/mapostolides/MetaFusion/RUNS
gene_bed=/Users/mapostolides/MetaFusion/reference_files/ens_known_genes.renamed.ENSG.bed
gene_info=/Users/mapostolides/MetaFusion/reference_files/Homo_sapiens.gene_info
genome_fasta=/Users/mapostolides/MetaFusion/reference_files/human_g1k_v37_decoy.fasta
recurrent_bedpe=/Users/mapostolides/MetaFusion/reference_files/blacklist_breakpoints.bedpe


date=July-30-2020
#DATASETS
sim45_sim52=1

# SIM45.SIM52.combined
if [ $sim45_sim52 -eq 1 ]; then
echo SIM45.SIM52
outdir=$test_dir/SIM45.SIM52.benchmark.$date.MetaFusion
cff=/Users/mapostolides/MetaFusion/test_data/cff/dream.sim45.sim52.cff
truth_fusions=/Users/mapostolides/MetaFusion/test_data/truth_sets/dream.sim45.sim52.truth_set.da

#sh MetaFusion.sh --outdir $outdir --cff $cff  --gene_bed $gene_bed  --gene_info $gene_info --truth_set $truth_fusions  --num_tools=2  --genome_fasta $genome_fasta --recurrent_bedpe $recurrent_bedpe
num_tools=2

sh MetaFusion.sh --outdir $outdir \
                 --cff $cff  \
                 --gene_bed $gene_bed \
                 --gene_info $gene_info \
                 --truth_set $truth_fusions \
                 --num_tools=2  \
                 --genome_fasta $genome_fasta \
                 --recurrent_bedpe $recurrent_bedpe \
                 --scripts $fusiontools
fi

