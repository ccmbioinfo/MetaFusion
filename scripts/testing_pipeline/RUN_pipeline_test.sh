#!/bin/bash
#NOTES
#$reann_test_dir/cff_iles/NOTES-merged.cff_RPRD2--LAMC2
#source /home/mapostolides/miniconda3/etc/profile.d/conda.sh
source /hpf/largeprojects/ccmbio/mapostolides/MODULES/miniconda3/etc/profile.d/conda.sh
conda activate metafusion

#CONSTANTS
#gene_info_file=/hpf/largeprojects/ccmbio/mapostolides/gene_fusion/pipeline/config_reference_files/Homo_sapiens.gene_info
info_file_dir=/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/rename_genes_test/info_file_subsets
reann_test_dir=/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/reann_cff_fusion_testing
test_dir=/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/testing_pipeline
outdir=$reann_test_dir/outdir-April-30-2020
#gene_bed_total=$reann_test_dir/ens_known_genes.renamed.bed
gene_bed_total=/hpf/largeprojects/ccmbio/mapostolides/gene_fusion/annotation/ens_known_genes.renamed.ENSG.bed
bed_dir=$reann_test_dir/bed_files
cff_dir=$reann_test_dir/cff_files
benchmark_toolkit=/hpf/largeprojects/ccmbio/mapostolides/MODULES/RUN_BENCHMARKING_TOOLKIT

run_pipeline (){
  #rename=1
  #annotate=1
  #merge=1
  #output_ANC_RT_SG=1
  RT_call_filter=1
  blck_filter=1
  ANC_filter=1
  benchmark=1
  outdir=$1
  mkdir -p $outdir
  cff=$2
  gene_bed=$3
  truth_fusions=$4
  num_tools=$5
  gene_info_file=/hpf/largeprojects/ccmbio/mapostolides/gene_fusion/pipeline/config_reference_files/Homo_sapiens.gene_info
  fusiontools=/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin
  genome_fasta=/hpf/largeprojects/ccmbio/mapostolides/gene_fusion/pipeline/config_reference_files/human_g1k_v37_decoy.fasta

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
    source /home/mapostolides/miniconda3/etc/profile.d/conda.sh
    sh $fusiontools/RUN_cluster_genes_breakpoints.sh $cff $outdir > $cluster
  fi
 
  #output ANC_RT_SG file
  if [ $output_ANC_RT_SG -eq 1 ]; then
    echo output ANC_RT_SG file
    python $fusiontools/output_ANC_RT_SG.py $cluster > $cluster.ANC_RT_SG 
  fi

  #ReadThrough Callerfilter
  if [ $RT_call_filter -eq 1 ]; then
    echo ReadThrough, callerfilter
    cat $cluster | grep ReadThrough > $outdir/$(basename $cluster).ReadThrough
    $fusiontools/callerfilter_num.py --cluster $cluster  --num_tools $num_tools | grep -v ReadThrough  > $outdir/$(basename $cluster).RT_filter.callerfilter.$num_tools
    #python $fusiontools/callerfilter_num.py --cluster $cluster  --num_tools $num_tools | grep -v ReadThrough  > $outdir/$(basename $cluster).RT_filter.callerfilter.$num_tools
    #cat $cluster | grep -v ReadThrough | awk '$8 ~ /.,./' > $outdir/$(basename $cluster).RT_filter.callerfilter2
  fi
  cluster_RT_call=$outdir/$(basename $cluster).RT_filter.callerfilter.$num_tools 

  # Blacklist Filter
  if [ $blck_filter -eq 1 ]; then
    echo blacklist filter
    blck_script_dir=/hpf/largeprojects/ccmbio/mapostolides/MODULES/FusionAnnotator/TEST_FusionAnnotator
  $blck_script_dir/blacklist_filter_recurrent_breakpoints.sh $cff $cluster_RT_call $outdir  > $outdir/$(basename $cluster).RT_filter.callerfilter.$num_tools.blck_filter
  fi
  cluster=$outdir/$(basename $cluster).RT_filter.callerfilter.$num_tools.blck_filter


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
    $benchmark_toolkit/benchmarking_cff_pertool.sh $outdir $truth_fusions $cff
}
benchmark_cff_renamed_reann_pertool (){
    cff=$1
    truth_fusions=$2
    outdir=$3
    mkdir -p $outdir
    # run benchmarking toolkit, no filters 
    $benchmark_toolkit/benchmarking_cff.renamed.reann_pertool.sh $outdir $truth_fusions $cff
}

run_fusionannotator () {
  outdir=$1
  cluster=$2
  echo Run FusionAnnotator and create cancer DB files
  /hpf/largeprojects/ccmbio/mapostolides/MODULES/RUN_BENCHMARKING_TOOLKIT/run_FusionAnnotator_cluster.sh $outdir $cluster
}

date=July-22-2020

#DATASETS
trusight=0
prostate=0
melanoma_cells=0
dipg=0
dipg_T=0
st_jude=0
brca_4=0
uhrr=0
beers_neg=0
sim_50=0
sim101=0
sim45_sim52=0
sim50_test=0

# ROB trusight
if [ $trusight -eq 1 ]; then
  cluster=/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/testing_pipeline/fusionOutput.txt
  outdir=$test_dir/ROB.TRUSIGHT.cancer_DB_hits.$date
  run_fusionannotator $outdir $cluster
fi

#SIM50 test cases
if [ $sim50_test -eq 1 ]; then
echo SIM50 test cases
outdir=$test_dir/SIM50.test_H2AJ.$date
truth_fusions=/hpf/largeprojects/ccmbio/mapostolides/MODULES/RUN_BENCHMARKING_TOOLKIT/renamed_truth_sets/sim50.truth_set.2500_fusions.dat
truth_fusions=/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/testing_pipeline/sim50.H2AJ--ADIG_truth.dat
cff=/hpf/largeprojects/ccmbio/mapostolides/validate_fusion/haas_2019_validation/output-SIM-April-17-2020/fusions/cff/merged.cff
#run_pipeline $outdir $cff $gene_bed $truth_fusions
cat $cff | sed 's/(.\+)//g' | sed 's/\//,/g' | grep arriba > $outdir/$(basename $cff).gene_names_cleaned.arriba
cff=$outdir/$(basename $cff).gene_names_cleaned
echo $cff
benchmark_pertool $cff $truth_fusions $outdir

fi

# Kumar prostate
if [ $prostate -eq 1 ]; then
echo Kumar prostate
gene_bed=$gene_bed_total
#outdir=$test_dir/kumar_prostate.benchmark.$date\.total_cluster
#outdir=$test_dir/kumar_prostate.benchmark.$date\.RTs_only
outdir=$test_dir/kumar_prostate.benchmark.with_chimerascan.$date
#outdir=$test_dir/kumar_prostate.benchmark.with_chimerascan.4callers.$date
mkdir -p $outdir
truth_fusions=/hpf/largeprojects/ccmbio/mapostolides/MODULES/RUN_BENCHMARKING_TOOLKIT/renamed_truth_sets/kumar.ReadThrough.truth_set.LNCaP_prostate_siCTCF.dat
# Run with CHIMERASCAN output
#cff=/hpf/largeprojects/ccmbio/mapostolides/MODULES/run_chimerascan/RUNS/prostate_cffs/8Callers.prostate.chimerascan_no_SRR1657556.cff
cff=/hpf/largeprojects/ccmbio/mapostolides/MODULES/run_chimerascan/RUNS/prostate_cffs/8Callers.prostate.chimerascan.cff
# Make all sample names the same
#cat $cff  | awk '{$8="LNCaP_prostate_siCTCF";print $0}'| sed 's/ /\t/g' > $outdir/$(basename $cff).prostate 
cff=$outdir/$(basename $cff).prostate
# Include only defuse ericscript integrate chimerascan
#cat $cff | grep 'defuse\|ericscript\|integrate\|chimerascan' > $outdir/$(basename $cff).4callers
#cff=$outdir/$(basename $cff).4callers

#run_pipeline $outdir $cff $gene_bed $truth_fusions
#Run benchmark on ANC_RT_SG .cluster
outdir=$test_dir/kumar_prostate.benchmark.chimerascan.ANC_RT_SG.$date
#outdir=$test_dir/kumar_prostate.benchmark.chimerascan.4callers.ANC_RT_SG.$date
mkdir $outdir
#cff_reann=/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/testing_pipeline/kumar_prostate.benchmark.with_chimerascan.July-2-2020/8Callers.prostate.chimerascan_no_SRR1657556.cff.prostate.renamed.reann
#cff_reann=/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/testing_pipeline/kumar_prostate.benchmark.with_chimerascan.July-7-2020/8Callers.prostate.chimerascan.cff.prostate.reformat.renamed.reann
#cff_reann=/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/testing_pipeline/kumar_prostate.benchmark.with_chimerascan.4callers.July-7-2020/8Callers.prostate.chimerascan.cff.prostate.4callers.reformat.renamed.reann
#cluster_ANC_RT_SG=/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/testing_pipeline/kumar_prostate.benchmark.with_chimerascan.July-2-2020/8Callers.prostate.chimerascan_no_SRR1657556.cff.prostate.renamed.reann.cluster.ANC_RT_SG
#cluster_ANC_RT_SG=/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/testing_pipeline/kumar_prostate.benchmark.with_chimerascan.July-7-2020/8Callers.prostate.chimerascan.cff.prostate.reformat.renamed.reann.cluster.ANC_RT_SG
#cluster_ANC_RT_SG=/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/testing_pipeline/kumar_prostate.benchmark.with_chimerascan.4callers.July-7-2020/8Callers.prostate.chimerascan.cff.prostate.4callers.reformat.renamed.reann.cluster.ANC_RT_SG

cff_reann=/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/testing_pipeline/kumar_prostate.benchmark.with_chimerascan.4callers.July-7-2020/8Callers.prostate.chimerascan.cff.prostate.4callers.reformat.renamed.reann
cluster_ANC_RT_SG=/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/testing_pipeline/kumar_prostate.benchmark.with_chimerascan.4callers.July-7-2020/8Callers.prostate.chimerascan.cff.prostate.4callers.reformat.renamed.reann.cluster.ANC_RT_SG
/hpf/largeprojects/ccmbio/mapostolides/MODULES/RUN_BENCHMARKING_TOOLKIT/benchmarking_cluster-GENAP.sh $outdir $truth_fusions $cff_reann $cluster_ANC_RT_SG true

#benchmark_pertool $cff $truth_fusions $outdir

#cluster_RT=/hpf/largeprojects/ccmbio/mapostolides/validate_fusion/kumar_2016_test/all-6-samples-pipeline-runs/output-June-17-2020/fusions/cff/merged.cff.renamed.reann.cluster.ReadThrough
#truth_fusions=/hpf/largeprojects/ccmbio/mapostolides/validate_fusion/kumar_2016_test/validation_testing/kumar.ReadThrough.truth_set.LNCaP_prostate_siCTCF.dat
# BENCHMARK ANC_RT_SG CLUSTER
#cff=/hpf/largeprojects/ccmbio/mapostolides/validate_fusion/kumar_2016_test/all-6-samples-pipeline-runs/output-June-17-2020/fusions/cff/merged.cff.renamed.reann
#cluster=/hpf/largeprojects/ccmbio/mapostolides/validate_fusion/kumar_2016_test/all-6-samples-pipeline-runs/output-June-17-2020/fusions/cff/merged.cff.renamed.reann.cluster
#cluster_filtered=/hpf/largeprojects/ccmbio/mapostolides/validate_fusion/kumar_2016_test/all-6-samples-pipeline-runs/output-June-17-2020/fusions/cff/merged.cff.renamed.reann.cluster.blck_filter.RT_filter.callerfilter2.ANC_filter
#python /hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/output_ANC_RT_SG.py $cluster > $outdir/$(basename $cluster).ANC_RT_SG
#cluster_ANC_RT_SG=$outdir/$(basename $cluster).ANC_RT_SG
#kumar et al paper don't provide samples corresponding to fusions. Replace sample name with common one
#cat $cluster_ANC_RT_SG  | awk '{$15="LNCaP_prostate_siCTCF";print $0}'| sed 's/ /\t/g' > $outdir/$(basename $cluster_ANC_RT_SG).prostate
#cluster_ANC_RT_SG=$outdir/$(basename $cluster_ANC_RT_SG).prostate
#cat $cluster_filtered  | awk '{$15="LNCaP_prostate_siCTCF";print $0}'| sed 's/ /\t/g' > $outdir/$(basename $cluster_filtered).prostate
#cluster_filtered=$outdir/$(basename $cluster_filtered).prostate
#/hpf/largeprojects/ccmbio/mapostolides/MODULES/RUN_BENCHMARKING_TOOLKIT/benchmarking_cluster-GENAP.sh $outdir $truth_fusions $cff $cluster_ANC_RT_SG true
#/hpf/largeprojects/ccmbio/mapostolides/MODULES/RUN_BENCHMARKING_TOOLKIT/benchmarking_cluster-GENAP.sh $outdir $truth_fusions $cff $cluster_filtered true

#kumar et al paper don't provide samples corresponding to fusions. Replace sample name with common one
#cat $cluster_RT  | awk '{$15="LNCaP_prostate_siCTCF";print $0}'| sed 's/ /\t/g' > $outdir/$(basename $cluster_RT).prostate
#cluster_RT=$outdir/$(basename $cluster_RT).prostate
#/hpf/largeprojects/ccmbio/mapostolides/MODULES/RUN_BENCHMARKING_TOOLKIT/benchmarking_cluster-GENAP.sh $outdir $truth_fusions $cff $cluster_RT true
#cff_orig=/hpf/largeprojects/ccmbio/mapostolides/validate_fusion/kumar_2016_test/all-6-samples-pipeline-runs/output-June-17-2020/fusions/cff/merged.cff
#cff_orig=$outdir/$(basename $cff).prostate
#cff_orig=$outdir/$(basename $cff).prostate.ReadThroughs
#cat $cff  | awk '{$8="LNCaP_prostate_siCTCF";print $0}'| sed 's/ /\t/g' > $cff_orig
#cat $cff  | grep ReadThrough | awk '{$8="LNCaP_prostate_siCTCF";print $0}'| sed 's/ /\t/g' > $cff_orig
#benchmark_pertool $cff_orig $truth_fusions $outdir
#END
fi

# Melanoma
if [ $melanoma_cells -eq 1 ]; then
echo Melanoma
gene_bed=$gene_bed_total
outdir=$test_dir/melanoma.benchmark.$date
truth_fusions=/hpf/largeprojects/ccmbio/mapostolides/validate_fusion/melanoma_validtaion/melanoma.truth_fusions.renamed.dat
cff=/hpf/largeprojects/ccmbio/mapostolides/validate_fusion/melanoma_validtaion/output-June-4-2020/fusions/cff/merged.cff

#caller #s 2-7
#run_pipeline $outdir $cff $gene_bed $truth_fusions 2
nums=$(echo 3 4 5 6 7)
for num in ${nums[@]};do
  cff=/hpf/largeprojects/ccmbio/mapostolides/validate_fusion/melanoma_validtaion/output-June-4-2020/fusions/cff/merged.cff
  run_pipeline $outdir $cff $gene_bed $truth_fusions $num
done
exit 0

#run_pipeline $outdir $cff $gene_bed $truth_fusions
#benchmark_pertool $cff $truth_fusions $outdir
#ADD cancer counts from Metacaller
meta_total=$(cat $outdir/merged.cff.renamed.reann.cluster.RT_filter.callerfilter2.blck_filter.ANC_filter | wc | awk '{print $1}')
meta_cancer=$(cat $outdir/merged.cff.renamed.reann.cluster.RT_filter.callerfilter2.blck_filter.ANC_filter.CANCER_FUSIONS | wc | awk '{print $1}')
echo -e MetaFuse\\t$meta_total\\t$meta_cancer >> $outdir/percaller.total.cancer_DB_hits
fi

# DIPG -- benchmark experimentally validated in RNA
if [ $dipg -eq 1 ]; then
echo DIPG
outdir=$test_dir/DIPG.RNA_validated.benchmark.$date
gene_bed=$gene_bed_total
truth_fusions=/hpf/largeprojects/ccmbio/mapostolides/MODULES/RUN_BENCHMARKING_TOOLKIT/renamed_truth_sets/DIPG.RNA_validated.renamed.NO_RT.dat
cff=/hpf/largeprojects/ccmbio/mapostolides/DIPG/run_DIPG_samples/DIPG_output_7CALLERS_April-17-2020/fusions/cff/merged.cff
run_pipeline $outdir $cff $gene_bed $truth_fusions
#Benchmark pertool
benchmark_pertool $cff $truth_fusions $outdir
# Benchmark cluster
cff=/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/testing_pipeline/DIPG.RNA_validated.benchmark.June-15-2020/merged.cff.renamed.reann
cluster=/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/testing_pipeline/DIPG.RNA_validated.benchmark.June-15-2020/merged.cff.renamed.reann.cluster.RT_filter.callerfilter2.blck_filter.ANC_filter
#/hpf/largeprojects/ccmbio/mapostolides/MODULES/RUN_BENCHMARKING_TOOLKIT/benchmarking_cluster-GENAP.sh $outdir $truth_fusions $cff $cluster true
#ADD cancer counts from Metacaller
meta_total=$(cat /hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/testing_pipeline/DIPG.RNA_validated.benchmark.June-15-2020/merged.cff.renamed.reann.cluster.RT_filter.callerfilter2.blck_filter.ANC_filter | wc | awk '{print $1}')
meta_cancer=$(cat /hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/testing_pipeline/DIPG.RNA_validated.benchmark.June-15-2020/merged.cff.renamed.reann.cluster.RT_filter.callerfilter2.blck_filter.ANC_filter.CANCER_FUSIONS | wc | awk '{print $1}')
echo -e MetaFuse\\t$meta_total\\t$meta_cancer >> $outdir/percaller.total.cancer_DB_hits
fi

# DIPG_T
if [ $dipg_T -eq 1 ]; then
echo DIPG_T
outdir=$test_dir/DIPG_T.RNA_validated.benchmark.$date
gene_bed=$gene_bed_total
truth_fusions=/hpf/largeprojects/ccmbio/mapostolides/MODULES/RUN_BENCHMARKING_TOOLKIT/renamed_truth_sets/DIPG.RNA_validated.renamed.NO_RT.dat
cff=/hpf/largeprojects/ccmbio/mapostolides/DIPG/run_DIPG_samples/DIPG_output_7CALLERS_April-17-2020/fusions/cff/merged.cff.T
run_pipeline $outdir $cff $gene_bed $truth_fusions
#Benchmark pertool
benchmark_pertool $cff $truth_fusions $outdir
# Benchmark cluster
#cff=/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/testing_pipeline/DIPG.RNA_validated.benchmark.June-15-2020/merged.cff.renamed.reann
#cluster=/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/testing_pipeline/DIPG.RNA_validated.benchmark.June-15-2020/merged.cff.renamed.reann.cluster.RT_filter.callerfilter2.blck_filter.ANC_filter
#/hpf/largeprojects/ccmbio/mapostolides/MODULES/RUN_BENCHMARKING_TOOLKIT/benchmarking_cluster-GENAP.sh $outdir $truth_fusions $cff $cluster true
#ADD cancer counts from Metacaller
#meta_total=$(cat /hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/testing_pipeline/DIPG.RNA_validated.benchmark.June-15-2020/merged.cff.renamed.reann.cluster.RT_filter.callerfilter2.blck_filter.ANC_filter | wc | awk '{print $1}')
#meta_cancer=$(cat /hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/testing_pipeline/DIPG.RNA_validated.benchmark.June-15-2020/merged.cff.renamed.reann.cluster.RT_filter.callerfilter2.blck_filter.ANC_filter.CANCER_FUSIONS | wc | awk '{print $1}')
#echo -e MetaFuse\\t$meta_total\\t$meta_cancer >> $outdir/percaller.total.cancer_DB_hits
fi

#STJUDE .cff fusions
if [ $st_jude -eq 1 ]; then
echo STJUDE
cff=/hpf/largeprojects/ccmbio/mapostolides/validate_fusion/stjude_validation/St-Jude-fusions.cff
outdir=$test_dir/STJUDE.annotate_cff_fusions.$date
gene_bed=$gene_bed_total
truth_fusions=/hpf/largeprojects/ccmbio/mapostolides/MODULES/RUN_BENCHMARKING_TOOLKIT/renamed_truth_sets/stjude.truth_set.no_NA.uniq.renamed.dat
#run_pipeline $outdir $cff $gene_bed $truth_fusions

#STJUDE
outdir=$test_dir/STJUDE.benchmark.June-12-2020.CodingFusion_truth.RT_filter.callerfilter2.blck_filter.RT_filter.callerfilter2.ANC_filter
#cluster=/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/testing_pipeline/STJUDE.benchmark.June-11-2020/merged.cff.renamed.reann.cluster.RT_filter.callerfilter2.blck_filter.RT_filter.callerfilter2.ANC_filter
#cluster=/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/testing_pipeline/STJUDE.benchmark.June-11-2020/merged.cff.renamed.reann.cluster
cluster=/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/testing_pipeline/STJUDE.benchmark.June-11-2020/merged.cff.renamed.reann.cluster.RT_filter.callerfilter2.blck_filter.RT_filter.callerfilter2.ANC_filter
cff=/hpf/largeprojects/ccmbio/mapostolides/validate_fusion/stjude_validation/stjude_output.April-23-2020/fusions/cff/merged.cff.renamed.reann
#truth_fusions=/hpf/largeprojects/ccmbio/mapostolides/MODULES/RUN_BENCHMARKING_TOOLKIT/renamed_truth_sets/stjude.truth_set.no_NA.uniq.renamed.dat
#truth_fusions=/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/testing_pipeline/STJUDE.annotate_cff_fusions.June-11-2020/stjude.truth_fusions.through_pipeline.dat
truth_fusions=/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/testing_pipeline/STJUDE.annotate_cff_fusions.June-11-2020/stjude.CodingFusion.truth_fusions.dat
#run_pipeline $outdir $cff $gene_bed $truth_fusions
#/hpf/largeprojects/ccmbio/mapostolides/MODULES/RUN_BENCHMARKING_TOOLKIT/benchmarking_cluster-GENAP.sh $outdir  $truth_fusions $cff $cluster
cff_orig=/hpf/largeprojects/ccmbio/mapostolides/validate_fusion/stjude_validation/stjude_output.April-23-2020/fusions/cff/merged.cff
#benchmark_pertool $cff_orig $truth_fusions $outdir
fi

#BT474.KPL4.MCF7.SKBR3
if [ $brca_4 -eq 1 ]; then
echo BT474.KPL4.MCF7.SKBR3
#outdir=$test_dir/BT474.KPL4.MCF7.SKBR3.benchmark.$date
#outdir=$test_dir/BT474.KPL4.MCF7.SKBR3.benchmark.June-18-2020
outdir=$test_dir/BT474.KPL4.MCF7.SKBR3.benchmark.$date
gene_bed=$gene_bed_total
truth_fusions=/hpf/largeprojects/ccmbio/mapostolides/MODULES/RUN_BENCHMARKING_TOOLKIT/renamed_truth_sets/BT474.KPL4.MCF7.SKBR3.truth_set.dat
cff=/hpf/largeprojects/ccmbio/mapostolides/validate_fusion/haas_2019_validation/output.BT474.KPL4.MCF7.SKBR3-April-9-2020/fusions/cff/merged.cff

#caller #s 2-7
nums=$(echo 2 3 4 5 6 7)
for num in ${nums[@]};do 
  cff=/hpf/largeprojects/ccmbio/mapostolides/validate_fusion/haas_2019_validation/output.BT474.KPL4.MCF7.SKBR3-April-9-2020/fusions/cff/merged.cff 
  run_pipeline $outdir $cff $gene_bed $truth_fusions $num
done
#run_pipeline $outdir $cff $gene_bed $truth_fusions 2
#run_pipeline $outdir $cff $gene_bed $truth_fusions 3
#run_pipeline $outdir $cff $gene_bed $truth_fusions 4
#run_pipeline $outdir $cff $gene_bed $truth_fusions 5
#run_pipeline $outdir $cff $gene_bed $truth_fusions 6
#run_pipeline $outdir $cff $gene_bed $truth_fusions 7

exit 0 

#benchmark_pertool $cff $truth_fusions $outdir
# benchmark renamed + annotated .cff, allows us to look at distribution of categories
cff_renamed_reann=/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/testing_pipeline/BT474.KPL4.MCF7.SKBR3.benchmark.June-18-2020/merged.cff.renamed.reann
#benchmark_cff_renamed_reann_pertool $cff_renamed_reann $truth_fusions $outdir

#Compare star_fusion FPs with MetaFusion FPs wrt cancer DB hits
outdir=$test_dir/BT474.KPL4.MCF7.SKBR3.star_fusion.MetaFusion.FP.cancer_hits-$date
cff_FP_star_fusion=$test_dir/BT474.KPL4.MCF7.SKBR3.benchmark.June-25-2020/merged.cff.renamed.reann.star_fusion.FP
cff_total=$test_dir/BT474.KPL4.MCF7.SKBR3.benchmark.June-29-2020/merged.cff.renamed.reann
cluster_FP=$test_dir/BT474.KPL4.MCF7.SKBR3.benchmark.June-29-2020/merged.cff.renamed.reann.cluster.RT_filter.callerfilter2.blck_filter.ANC_filter.FP
$benchmark_toolkit/benchmarking_cluster-GENAP.sh $outdir $truth_fusions $cff_total $cluster_FP true
$benchmark_toolkit/benchmarking_cff_pertool.sh $outdir $truth_fusions $cff_FP_star_fusion
fi

#UHRR
if [ $uhrr -eq 1 ]; then
echo UHRR
#cat merged.cff.renamed.reann.cluster.blck_filter.RT_filter.callerfilter2.ANC_filter | awk '{$15="UHR";print $0}'| sed 's/ /\t/g'  > merged.cff.renamed.reann.cluster.blck_filter.RT_filter.callerfilter2.ANC_filter.UHRR; mv merged.cff.renamed.reann.cluster.blck_filter.RT_filter.callerfilter2.ANC_filter.UHRR merged.cff.renamed.reann.cluster.blck_filter.RT_filter.callerfilter2.ANC_filter
outdir=$test_dir/UHRR.benchmark.$date
mkdir $outdir
gene_bed=$gene_bed_total
cff=/hpf/largeprojects/ccmbio/mapostolides/fusion_pipeline_runs/UHRR/output_Feb5_2020/fusions/cff/merged.cff
#rename samples to UHR
cat $cff  | awk '{$8="UHR";print $0}'| sed 's/ /\t/g' > $outdir/$(basename $cff).UHR
cff=$outdir/$(basename $cff).UHR
#cluster=/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/testing_pipeline/UHRR.benchmark.June-9-2020/merged.cff.renamed.reann.cluster.blck_filter.RT_filter.callerfilter2.ANC_filter
truth_fusions=/hpf/largeprojects/ccmbio/mapostolides/MODULES/RUN_BENCHMARKING_TOOLKIT/renamed_truth_sets/UHR.truth_fusions.renamed.dat
#sh /hpf/largeprojects/ccmbio/mapostolides/MODULES/RUN_BENCHMARKING_TOOLKIT/benchmarking_cluster-GENAP.sh $outdir $truth_fusions $cff $cluster true 
run_pipeline $outdir $cff $gene_bed $truth_fusions
#cff_ann=$outdir/$cff.reann
benchmark_pertool $cff $truth_fusions $outdir
#ADD cancer counts from Metacaller
meta_total=$(cat $outdir/merged.cff.UHR.renamed.reann.cluster.RT_filter.callerfilter2.blck_filter.ANC_filter | wc | awk '{print $1}')
meta_cancer=$(cat $outdir/merged.cff.UHR.renamed.reann.cluster.RT_filter.callerfilter2.blck_filter.ANC_filter.CANCER_FUSIONS | wc | awk '{print $1}')
echo -e MetaFuse\\t$meta_total\\t$meta_cancer >> $outdir/percaller.total.cancer_DB_hits
fi

#NEGATIVE CONTROL BEERS
if [ $beers_neg -eq 1 ]; then
echo NEGATIVE CONTROL BEERS
outdir=$test_dir/NEG_CONTROL_BEERS.benchmark.$date
gene_bed=$gene_bed_total
cff=/hpf/largeprojects/ccmbio/mapostolides/validate_fusion/montreal_benchmark_validation/negative_control_beers/output-7CALLERS-April-17-2020/fusions/cff/merged.cff
truth_fusions=/hpf/largeprojects/ccmbio/mapostolides/MODULES/RUN_BENCHMARKING_TOOLKIT/renamed_truth_sets/DIPG.RNA_validated.renamed.dat
#run_pipeline $outdir $cff $gene_bed $truth_fusions
fi

# SIM50 2500 fusions files:
if [ $sim_50 -eq 1 ]; then
echo SIM50
outdir=$test_dir/SIM50.2500_TP.benchmark.$date
#outdir=$test_dir/SIM50.benchmark_pertool.benchmark_HGNC_updated.$date
gene_bed=$gene_bed_total
truth_fusions=/hpf/largeprojects/ccmbio/mapostolides/MODULES/RUN_BENCHMARKING_TOOLKIT/renamed_truth_sets/sim50.truth_set.2500_fusions.dat
#truth_fusions=/hpf/largeprojects/ccmbio/mapostolides/validate_fusion/haas_2019_validation/haas_2019_simulated_dataset_FILES/sim_50.truth_set.dat
cff=/hpf/largeprojects/ccmbio/mapostolides/validate_fusion/haas_2019_validation/output-SIM-April-17-2020/fusions/cff/merged.cff

#caller #s 2-7
#run_pipeline $outdir $cff $gene_bed $truth_fusions 2
nums=$(echo 3 4 5 6 7)
for num in ${nums[@]};do
  cff=/hpf/largeprojects/ccmbio/mapostolides/validate_fusion/haas_2019_validation/output-SIM-April-17-2020/fusions/cff/merged.cff
  run_pipeline $outdir $cff $gene_bed $truth_fusions $num
done
exit 0

run_pipeline $outdir $cff $gene_bed $truth_fusions
echo $cff
#benchmark_pertool $cff $truth_fusions $outdir
cff_renamed_reann=$outdir/merged.cff.renamed.reann
benchmark_cff_renamed_reann_pertool $cff_renamed_reann $truth_fusions $outdir
fi

#SIM101 2500 fusions files, same truth set as SIM50
if [ $sim101 -eq 1 ]; then
echo SIM101
outdir=$test_dir/SIM101.2500_TP.benchmark.$date
truth_fusions=/hpf/largeprojects/ccmbio/mapostolides/MODULES/RUN_BENCHMARKING_TOOLKIT/renamed_truth_sets/sim_101.truth_set.2500_fusions.dat
gene_bed=$gene_bed_total
cff=/hpf/largeprojects/ccmbio/mapostolides/validate_fusion/haas_2019_validation/output-SIM_101-April-23-2020/fusions/cff/merged.cff

#caller #s 2-7
#run_pipeline $outdir $cff $gene_bed $truth_fusions 2
nums=$(echo 3 4 5 6 7)
for num in ${nums[@]};do
  cff=/hpf/largeprojects/ccmbio/mapostolides/validate_fusion/haas_2019_validation/output-SIM_101-April-23-2020/fusions/cff/merged.cff
  run_pipeline $outdir $cff $gene_bed $truth_fusions $num
done
exit 0

run_pipeline $outdir $cff $gene_bed $truth_fusions
#benchmark_pertool $cff $truth_fusions $outdir
cff_renamed_reann=$outdir/merged.cff.renamed.reann
benchmark_cff_renamed_reann_pertool $cff_renamed_reann $truth_fusions $outdir
fi

# SIM45.SIM52.combined
if [ $sim45_sim52 -eq 1 ]; then
echo SIM45.SIM52
gene_bed=$gene_bed_total
outdir=$test_dir/SIM45.SIM52.benchmark.$date
cff=/hpf/largeprojects/ccmbio/mapostolides/validate_fusion/montreal_benchmark_validation/sim45.sim52.combined/output-7CALLERS-April-22-2020/fusions/cff/merged.cff
truth_fusions=/hpf/largeprojects/ccmbio/mapostolides/MODULES/RUN_BENCHMARKING_TOOLKIT/renamed_truth_sets/sim45.sim52.truth_set.dat

#caller #s 2-7
#run_pipeline $outdir $cff $gene_bed $truth_fusions 2
nums=$(echo 3 4 5 6 7)
for num in ${nums[@]};do
  cff=/hpf/largeprojects/ccmbio/mapostolides/validate_fusion/montreal_benchmark_validation/sim45.sim52.combined/output-7CALLERS-April-22-2020/fusions/cff/merged.cff
  run_pipeline $outdir $cff $gene_bed $truth_fusions $num
done
exit 0


run_pipeline $outdir $cff $gene_bed $truth_fusions
benchmark_pertool $cff $truth_fusions $outdir
fi

exit 0
# TEST FUSIONS BELOW THIS POINT
outdir=$test_dir/outdir.TEST_FUSIONS.May-20-2020
#BCL3--CTB-171A8.1
#/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/reann_cff_fusion_testing/cff_files/BCL3--CTB-171A8.1.cff
cff=$cff_dir/BCL3--CTB-171A8.1.cff
gene_bed=$bed_dir/BCL3--CTB-171A8.1.bed
info_file=$info_file_dir/Homo_sapiens.gene_info.BCL3
#run_pipeline $outdir $cff $gene_bed


echo ARRIBA
cff=/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/rename_genes_test/arriba.weird_brackets.cff
#python $fusiontools/rename_cff_file_genes-GENAP.py $cff $gene_info_file
run_pipeline $outdir $cff $gene_bed
echo INTEGRATE
cff=/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/rename_genes_test/integrate_slash.delimeter.cff
run_pipeline $outdir $cff $gene_bed
#python $fusiontools/rename_cff_file_genes-GENAP.py $cff $gene_info_file

# TEST FUSIONS
#HES1--IL29
cff=$reann_test_dir/HES1--IL29/HES1--IL29.cff
gene_bed=$reann_test_dir/ens_known_genes_HES1_IL29.bed
run_pipeline $outdir $cff $gene_bed
#SAA1--ATP1A1
cff=$reann_test_dir/SAA1--ATP1A1/SAA1--ATP1A1.cff
gene_bed=$reann_test_dir/ens_known_genes_SAA1--ATP1A1.bed
run_pipeline $outdir $cff $gene_bed
#BCL3--CTB-171A8.1
cff=$cff_dir/BCL3--CTB-171A8.1.cff
gene_bed=$bed_dir/BCL3--CTB-171A8.1.bed
run_pipeline $outdir $cff $gene_bed
#ANKIB1--AKAP9
cff=$reann_test_dir/cff_files/merged.cff_AKAP9--ANKIB1
gene_bed=$reann_test_dir/ens_known_genes.AKAP9--ANKIB1.bed
run_pipeline $outdir $cff $gene_bed
#FOSB--AADACL2
cff=$reann_test_dir/FOSB--AADACL2/FOSB--AADACL2.cff
gene_bed=$bed_dir/ens_known_genes.MIR548H2.FOSB.AADACL2.bed
run_pipeline $outdir $cff $gene_bed


#SAA1--ATP1A1
outdir=$reann_test_dir/SAA1--ATP1A1
cff=$reann_test_dir/SAA1--ATP1A1/SAA1--ATP1A1.cff
#run_pipeline $outdir $cff $gene_bed_total

##
cff=$reann_test_dir/cff_files/merged.cff_BCL3 

#BCL3--CTB-171A8.1
cff=$cff_dir/BCL3--CTB-171A8.1.cff
gene_bed=$bed_dir/BCL3--CTB-171A8.1.bed
gene_bed=$gene_bed_total
#ANKIB1--AKAP9
#gene_bed=$reann_test_dir/ens_known_genes.ANKIB1--AKAP9.bed
#cff=$reann_test_dir/cff_files/merged.cff_AKAP9--ANKIB1
#run_pipeline $outdir $cff $gene_bed

#FOSB--AADACL2
cff=$reann_test_dir/FOSB--AADACL2/FOSB--AADACL2.cff
gene_bed=$bed_dir/ens_known_genes.MIR548H2.FOSB.AADACL2.bed
#1 entry, failing due to breakpt being 1 bp off:
#cff=$reann_test_dir/FOSB--AADACL2/FOSB--AADACL2.cff.test
#run_pipeline $outdir $cff $gene_bed

sh /hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/testing_pipeline/testing.sh -a false -b false -f /hpf/largeprojects/ccmbio/mapostolides/MODULES/RUN_BENCHMARKING_TOOLKIT/run_FusionAnnotator_cluster.sh -v false

