#!/bin/bash
#NOTES
#source /hpf/largeprojects/ccmbio/mapostolides/MODULES/miniconda3/etc/profile.d/conda.sh
#conda activate metafusion

fusiontools=/MetaFusion/scripts
#REFERENCE FILES FILES
test_dir=/MetaFusion/RUNS
gene_bed=/MetaFusion/reference_files/ens_known_genes.renamed.ENSG.bed
gene_info=/MetaFusion/reference_files/Homo_sapiens.gene_info
genome_fasta=/MetaFusion/reference_files/human_g1k_v37_decoy.fasta
recurrent_bedpe=/MetaFusion/reference_files/blacklist_breakpoints.bedpe


date=Aug-6-2020.WITH_SEQ
#DATASETS
#sim45_sim52=1
brca_4=1
beers_neg=1
sim_50=1
sim101=1

# SIM45.SIM52.combined
if [ $sim45_sim52 -eq 1 ]; then
echo SIM45.SIM52
outdir=$test_dir/SIM45.SIM52.benchmark.$date.MetaFusion
cff=/MetaFusion/test_data/cff/dream.sim45.sim52.cff
truth_fusions=/MetaFusion/test_data/truth_sets/dream.sim45.sim52.truth_set.dat

#sh MetaFusion.sh --outdir $outdir --cff $cff  --gene_bed $gene_bed  --gene_info $gene_info --truth_set $truth_fusions  --num_tools=2  --genome_fasta $genome_fasta --recurrent_bedpe $recurrent_bedpe
bash MetaFusion.sh --outdir $outdir \
                 --cff $cff  \
                 --gene_bed $gene_bed \
                 --gene_info $gene_info \
                 --truth_set $truth_fusions \
                 --genome_fasta $genome_fasta \
                 --num_tools=2  \
                 --recurrent_bedpe $recurrent_bedpe \
                 --scripts $fusiontools
fi

#BT474.KPL4.MCF7.SKBR3
if [ $brca_4 -eq 1 ]; then
echo BT474.KPL4.MCF7.SKBR3
outdir=$test_dir/BT474.KPL4.MCF7.SKBR3.benchmark.$date
truth_fusions=/MetaFusion/test_data/truth_sets/BRCA.truth_set.dat
cff=/MetaFusion/test_data/cff/BRCA.cff

bash MetaFusion.sh --outdir $outdir \
                 --cff $cff  \
                 --gene_bed $gene_bed \
                 --gene_info $gene_info \
                 --genome_fasta $genome_fasta \
                 --truth_set $truth_fusions \
                 --num_tools=2  \
                 --recurrent_bedpe $recurrent_bedpe \
                 --scripts $fusiontools
fi

#NEGATIVE CONTROL BEERS
if [ $beers_neg -eq 1 ]; then
echo NEGATIVE CONTROL BEERS
outdir=$test_dir/NEG_CONTROL_BEERS.benchmark.$date
cff=/MetaFusion/test_data/cff/beers_neg.cff 
truth_fusions=/MetaFusion/test_data/truth_sets/BRCA.truth_set.dat

bash MetaFusion.sh --outdir $outdir \
                 --cff $cff  \
                 --gene_bed $gene_bed \
                 --gene_info $gene_info \
                 --genome_fasta $genome_fasta \
                 --truth_set $truth_fusions \
                 --num_tools=2  \
                 --recurrent_bedpe $recurrent_bedpe \
                 --scripts $fusiontools
fi

# SIM50 2500 fusions files:
if [ $sim_50 -eq 1 ]; then
echo SIM50
outdir=$test_dir/SIM50.2500_TP.benchmark.$date
cff=/MetaFusion/test_data/cff/sim50.cff
truth_fusions=/MetaFusion/test_data/truth_sets/sim50.truth_set.dat

bash MetaFusion.sh --outdir $outdir \
                 --cff $cff  \
                 --gene_bed $gene_bed \
                 --gene_info $gene_info \
                 --genome_fasta $genome_fasta \
                 --truth_set $truth_fusions \
                 --num_tools=2  \
                 --recurrent_bedpe $recurrent_bedpe \
                 --scripts $fusiontools

fi

#SIM101 2500 fusions files, same truth set as SIM50
if [ $sim101 -eq 1 ]; then
echo SIM101
outdir=$test_dir/SIM101.2500_TP.benchmark.$date
cff=/MetaFusion/test_data/cff/sim101.cff
truth_fusions=/MetaFusion/test_data/truth_sets/sim101.truth_set.dat

bash MetaFusion.sh --outdir $outdir \
                 --cff $cff  \
                 --gene_bed $gene_bed \
                 --gene_info $gene_info \
                 --genome_fasta $genome_fasta \
                 --truth_set $truth_fusions \
                 --num_tools=2  \
                 --recurrent_bedpe $recurrent_bedpe \
                 --scripts $fusiontools

fi



