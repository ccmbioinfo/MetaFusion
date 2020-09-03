#!/bin/bash

#Change date to current date. Can also add tag to this string for multiple runs
date=Sept-1-2020

fusiontools=/MetaFusion/scripts
#REFERENCE FILES FILES
runs_dir=/MetaFusion/EXON_RUNS
mkdir $runs_dir
gene_bed=/MetaFusion/reference_files/ens_known_genes.renamed.ENSG.bed
gene_info=/MetaFusion/reference_files/Homo_sapiens.gene_info
genome_fasta=/MetaFusion/reference_files/human_g1k_v37_decoy.fasta
recurrent_bedpe=/MetaFusion/reference_files/blacklist_breakpoints.bedpe

#DATASETS
#notch_nup=1
#stx_rae=1
#pvt_pdg=1
brca_4=1

extract_exon_pipeline (){
  outdir=$1
  mkdir -p $outdir
  cff=$2
  gene_bed=$3
  genome_fasta=$4
  
  fusiontools=/MetaFusion/scripts 
  #REFORMAT
  cat $cff | awk '$1 ~ /[0-9XY]/ && $4 ~ /[0-9XY]/ ' |  awk 'BEGIN{FS=OFS="\t"} $3 !~ /^[-+]$/{$3="NA"} 1' | awk 'BEGIN{FS=OFS="\t"} $6 !~ /^[-+]$/{$6="NA"} 1'   > $outdir/$(basename $cff).reformat
  cff=$outdir/$(basename $cff).reformat

  echo Annotate cff, no extraction of sequence surrounding breakpoint
  python reann_cff_fusion.py --cff $cff --gene_bed $gene_bed > $outdir/$(basename $cff).reann.NO_SEQ
  cff=$outdir/$(basename $cff).reann.NO_SEQ

  echo Add adjacent exons to cff
  python $fusiontools/extract_closest_exons.py $cff $gene_bed $genome_fasta  > $outdir/$(basename $cff).exons
  cff=$outdir/$(basename $cff).exons

  echo Merge cff by exon
  cluster=$outdir/$(basename $cff).cluster
  bash RUN_cluster_exons.sh $cff $outdir $fusiontools > $cluster
  #ReadThrough Callerfilter

  echo ReadThrough, callerfilter
  num_tools=2
  python callerfilter_num.py --cluster $cluster  --num_tools $num_tools | grep -v ReadThrough  > $outdir/$(basename $cluster).RT_filter.callerfilter.$num_tools
  cluster_RT_call=$outdir/$(basename $cluster).RT_filter.callerfilter.$num_tools

  # Blocklist Filter
  echo blocklist filter
  bash blocklist_filter_recurrent_breakpoints.sh $cff $cluster_RT_call $outdir $recurrent_bedpe > $outdir/$(basename $cluster).RT_filter.callerfilter.$num_tools.blck_filter
  cluster=$outdir/$(basename $cluster).RT_filter.callerfilter.$num_tools.blck_filter

  echo ANC adjacent noncoding filter
  python filter_adjacent_noncoding.py $cluster > $outdir/$(basename $cluster).ANC_filter
  cluster=$outdir/$(basename $cluster).ANC_filter

  echo Rank and generate final.cluster
  python rank_cluster_file.py $cluster > $outdir/final.cluster
  cluster=$outdir/final.cluster

  #echo Running benchmarking
  #echo $outdir
  #truth_set=/MetaFusion/test_data/truth_sets/BRCA.truth_set.dat
  #echo $truth_set
  #echo $cluster
  #echo $fusiontools
  #bash benchmarking_cluster.MetaFusion.sh $outdir $truth_set $cluster $fusiontools

}

if [ $notch_nup -eq 1 ]; then
#NOTCH1--NUP214
echo NOTCH1--NUP214
cff=/MetaFusion/test_data/cff_test/NOTCH1--NUP214.cff
outdir=$runs_dir/outdir-$date
mkdir $outdir
gene_bed=/MetaFusion/test_data/cff_test/ens_known_genes.NOTCH1--NUP214.bed
extract_exon_pipeline $outdir $cff $gene_bed $genome_fasta 
#python $fusiontools/reann_cff_fusion_exon_test2.py $cff $gene_bed $genome_fasta > $outdir/$(basename $cff).exons
fi

if [ $stx_rae -eq 1 ]; then
#STX16--RAE1
cff=/MetaFusion/test_data/cff/STX16--RAE1.ALL.cff
#cff=/MetaFusion/EXON_RUNS/outdir-Aug-27-2020/STX16--RAE1.ALL.cff.reann.NO_SEQ.55929088.arriba.1entry
outdir=$runs_dir/STX16--RAE1-$date
mkdir $outdir
gene_bed=/MetaFusion/test_data/cff_test/ens_known_genes.STX16--RAE1.renamed.bed
extract_exon_pipeline $outdir $cff $gene_bed $genome_fasta 
fi

if [ $pvt_pdg -eq 1 ]; then
#PVT1--PDGFB
cff=/MetaFusion/test_data/cff_test/PVT1--PDGFB.cff
outdir=$runs_dir/outdir-$date
mkdir $outdir
gene_bed=/MetaFusion/test_data/cff_test/ens_known_genes.renamed.PVT1--PDGFB.ENSG.bed
extract_exon_pipeline $outdir $cff $gene_bed $genome_fasta 
fi

#BT474.KPL4.MCF7.SKBR3
if [ $brca_4 -eq 1 ]; then
echo BT474.KPL4.MCF7.SKBR3
outdir=$runs_dir/BT474.KPL4.MCF7.SKBR3.exon_annotation.$date
echo generating output in $outdir
cff=/MetaFusion/test_data/cff/BRCA.cff
extract_exon_pipeline $outdir $cff $gene_bed $genome_fasta
fi

exit 0
#SCRAP

# SIM45.SIM52.combined
if [ $sim45_sim52 -eq 1 ]; then
echo SIM45.SIM52
outdir=$runs_dir/SIM45.SIM52.benchmark.$date
echo generating output in $outdir
mkdir $outdir
cff=/MetaFusion/test_data/cff/dream.sim45.sim52.cff

fi

#NEGATIVE CONTROL BEERS
if [ $beers_neg -eq 1 ]; then
echo NEGATIVE CONTROL BEERS
outdir=$runs_dir/BEERS.$date
echo generating output in $outdir
cff=/MetaFusion/test_data/cff/beers_neg.cff 
truth_fusions=/MetaFusion/test_data/truth_sets/BRCA.truth_set.dat

fi

# SIM50 2500 fusions files:
if [ $sim_50 -eq 1 ]; then
echo SIM50
outdir=$runs_dir/SIM50.$date
cff=/MetaFusion/test_data/cff/sim50.cff
fi

#SIM101 2500 fusions files, same truth set as SIM50
if [ $sim101 -eq 1 ]; then
echo SIM101
outdir=$runs_dir/SIM101.$date
echo generating output in $outdir
cff=/MetaFusion/test_data/cff/sim101.cff
fi

# STX16--RAE1 from BT474.KPL4.MCF7.SKBR3 dataset
if [ $stx16_rae1 -eq 1 ]; then
echo STX16--RAE1
outdir=$runs_dir/STX16--RAE1.$date
echo generating output in $outdir
cff=/MetaFusion/test_data/cff/STX16--RAE1.figure_subset.cff
fi

# Melanoma and CML
if [ $melanoma -eq 1 ]; then
echo MELANOMA and CML 
outdir=$runs_dir/melanoma.CML.$date
cff=/MetaFusion/test_data/cff/melanoma.cff
fi

