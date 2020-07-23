#!/bin/bash
#NOTES
#$reann_test_dir/cff_files/NOTES-merged.cff_RPRD2--LAMC2
#source /home/mapostolides/miniconda3/etc/profile.d/conda.sh
#source /hpf/largeprojects/ccmbio/mapostolides/MODULES/miniconda3/etc/profile.d/conda.sh
#conda activate metafusion

#CONSTANTS
run_dir=/Users/mapostolides/MetaFusion/RUNS
#gene_bed_total=$reann_test_dir/ens_known_genes.renamed.bed
gene_bed_total=/Users/mapostolides/MetaFusion/reference_files/ens_known_genes.renamed.ENSG.bed

run_pipeline (){
  rename=1
  annotate=1
  merge=1
  output_ANC_RT_SG=1
  RT_call_filter=1
  blck_filter=1
  ANC_filter=1
  #benchmark=0
  outdir=$1
  mkdir -p $outdir
  cff=$2
  gene_bed=$3
  gene_info_file=/Users/mapostolides/MetaFusion/reference_files/Homo_sapiens.gene_info
  truth_fusions=$4
  fusiontools=/Users/mapostolides/MetaFusion/mugqic_tools/python-tools/fusiontools/0.1.0/bin
  genome_fasta=/Users/mapostolides/MetaFusion/reference_files/human_g1k_v37_decoy.fasta
  recurrent_bedpe=/Users/mapostolides/MetaFusion/reference_files/blacklist_breakpoints.bedpe

  #Check CFF file format:
  #Remove entries with nonconformming chromosome name
  cat $cff | awk '$1 ~ /[0-9XY]/ && $4 ~ /[0-9XY]/ ' > $outdir/$(basename $cff).reformat 
  cff=$outdir/$(basename $cff).reformat
  #NEED TO INSERT +/-/NA for strand, make NA if other

  #Rename cff
  if [ $rename -eq 1 ]; then
    echo Rename cff
    python $fusiontools/rename_cff_file_genes-GENAP.py $cff $gene_info_file > $outdir/$(basename $cff).renamed
  fi
  cff=$outdir/$(basename $cff).renamed

  #Annotate cff
  if [ $annotate -eq 1 ]; then
    echo Annotate cff
    python $fusiontools/reann_cff_fusion.py $cff $gene_bed $genome_fasta > $outdir/$(basename $cff).reann 
  fi
  cff=$outdir/$(basename $cff).reann

  #Merge
  cluster=$outdir/$(basename $cff).cluster
  if [ $merge -eq 1 ]; then
    echo Merge cff
    #source /home/mapostolides/miniconda3/etc/profile.d/conda.sh
    #sh $fusiontools/RUN_cluster_genes_breakpoints.sh $cff $outdir > $cluster
    sh $fusiontools/RUN_cluster_genes_breakpoints.local.sh $cff $outdir $fusiontools > $cluster
  fi
 
  #output ANC_RT_SG file
  if [ $output_ANC_RT_SG -eq 1 ]; then
    echo output ANC_RT_SG file
    python $fusiontools/output_ANC_RT_SG.py $cluster > $cluster.ANC_RT_SG 
  fi

  #ReadThrough Callerfilter2
  if [ $RT_call_filter -eq 1 ]; then
    echo ReadThrough, callerfilter2
    cat $cluster | grep ReadThrough > $outdir/$(basename $cluster).ReadThrough
    cat $cluster | grep -v ReadThrough | awk '$8 ~ /.,./' > $outdir/$(basename $cluster).RT_filter.callerfilter2
    #cat $cluster | grep -v ReadThrough  > $outdir/$(basename $cluster).RT_filter.callerfilter2
  fi
  cluster_RT_call=$outdir/$(basename $cluster).RT_filter.callerfilter2 

  # Blacklist Filter
  if [ $blck_filter -eq 1 ]; then
    echo blacklist filter
    $fusiontools/blacklist_filter_recurrent_breakpoints.local.sh $cff $cluster_RT_call $outdir $recurrent_bedpe  > $outdir/$(basename $cluster).RT_filter.callerfilter2.blck_filter
  fi
  cluster=$outdir/$(basename $cluster).RT_filter.callerfilter2.blck_filter

  # Adjacent Noncoding filter 
  if [ $ANC_filter -eq 1 ]; then
    echo ANC adjacent noncoding filter
    $fusiontools/filter_adjacent_noncoding.py $cluster > $outdir/$(basename $cluster).ANC_filter  
  fi
  cluster=$outdir/$(basename $cluster).ANC_filter

  #Benchmark
  if [ $benchmark -eq 1 ]; then
     echo benchmark
    /hpf/largeprojects/ccmbio/mapostolides/MODULES/RUN_BENCHMARKING_TOOLKIT/benchmarking_cluster-GENAP.sh $outdir $truth_fusions $cff $cluster true 
  fi
}
benchmark_pertool (){
    cff=$1
    truth_fusions=$2
    outdir=$3
    mkdir -p $outdir
    # run benchmarking toolkit, no filters 
    /hpf/largeprojects/ccmbio/mapostolides/MODULES/RUN_BENCHMARKING_TOOLKIT/benchmarking_cff_pertool.sh $outdir $truth_fusions $cff
}
benchmark_cff_renamed_reann_pertool (){
    cff=$1
    truth_fusions=$2
    outdir=$3
    mkdir -p $outdir
    # run benchmarking toolkit, no filters 
    /hpf/largeprojects/ccmbio/mapostolides/MODULES/RUN_BENCHMARKING_TOOLKIT/benchmarking_cff.renamed.reann_pertool.sh $outdir $truth_fusions $cff
}

date=July-10-2020
#DATASETS
echo SIM45.SIM52
gene_bed=$gene_bed_total
outdir=$run_dir/SIM45.SIM52.$date
cff=/Users/mapostolides/MetaFusion/RUNS/cff_files/sim45.sim52.merged.cff
truth_fusions=/hpf/largeprojects/ccmbio/mapostolides/MODULES/RUN_BENCHMARKING_TOOLKIT/renamed_truth_sets/sim45.sim52.truth_set.dat
run_pipeline $outdir $cff $gene_bed $truth_fusions
#benchmark_pertool $cff $truth_fusions $outdir


