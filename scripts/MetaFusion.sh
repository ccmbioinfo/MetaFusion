#!/bin/bash
#outdir=$1
#mkdir -p $outdir
#cff=$2
#gene_bed=$3
#truth_fusions=$4
#num_tools=$5
#gene_info=$6
#genome_fasta=$7

#STEPS
rename=1
annotate=1
merge=1
output_ANC_RT_SG=1
RT_call_filter=1
blck_filter=1
ANC_filter=1
rank=1
benchmark=1


# Loop through arguments and process them
#for arg in "$@"
while test $# -gt 0;do
    case $1 in
        -n=*|--num_tools=*)
        num_tools="${1#*=}"
        shift
        ;;
        --outdir)
        outdir="$2"
        shift 2
        ;;
        --cff)
        cff="$2"
        shift 2
        ;;
        --gene_bed)
        gene_bed="$2"
        shift 2
        ;;
        --gene_info)
        gene_info="$2"
        shift 2
        ;;
        --genome_fasta)
        genome_fasta="$2"
        shift 2
        ;;
        --truth_set)
        truth_set="$2"
        shift 2
        ;;
        --recurrent_bedpe)
        recurrent_bedpe="$2"
        shift 2
        ;;
        --scripts)
        fusiontools="$2"
        shift 2
        ;;
        --fusion_annotator)
        FA=1
        shift
        ;;
        *)
        #OTHER_ARGUMENTS+=("$1")
        shift # Remove generic argument from processing
        ;;
    esac
done
#echo num_tools $num_tools 
#echo outdir $outdir 
#echo cff  $cff 
#echo gene_bed $gene_bed 
#echo gene_info $gene_info 
#echo genome_fasta $genome_fasta 
#echo truth_set $truth_set

#METAFUSION WORKFLOW
mkdir $outdir

#Check CFF file format:
#Remove entries with nonconformming chromosome name
#Remove "." from strand field and replace with "NA"
cat $cff | awk '$1 ~ /[0-9XY]/ && $4 ~ /[0-9XY]/ ' |  awk 'BEGIN{FS=OFS="\t"} $3 !~ /^[-+]$/{$3="NA"} 1' | awk 'BEGIN{FS=OFS="\t"} $6 !~ /^[-+]$/{$6="NA"} 1'   > $outdir/$(basename $cff).reformat 
#cat $cff | awk '$1 ~ /[0-9XY]/ && $4 ~ /[0-9XY]/ ' |  awk 'FS=OFS="\t"{if ($3==".") {$3="NA"}; print}'|  awk 'FS=OFS="\t"{if ($6==".") {$6="NA"}; print}'  > $outdir/$(basename $cff).reformat 
#cat $cff | awk '$1 ~ /[0-9XY]/ && $4 ~ /[0-9XY]/ ' | awk 'FS=OFS="\t"{ if ($3 !="+" || $3 != "-" ); $3="NA"; print }' | awk 'FS=OFS="\t"{ if ($6 !="+" || $6 != "-" ); $6="NA"; print }' > $outdir/$(basename $cff).reformat 
cff=$outdir/$(basename $cff).reformat
#NEED TO INSERT +/-/NA for strand, make NA if other


#Rename cff
if [ $rename -eq 1 ]; then
  echo Rename cff
  python rename_cff_file_genes-GENAP.py $cff $gene_info > $outdir/$(basename $cff).renamed
fi
cff=$outdir/$(basename $cff).renamed

#Annotate cff
if [ $annotate -eq 1 ]; then
  if [ $genome_fasta ]; then 
    echo Annotate cff, extract sequence surrounding breakpoint
    python reann_cff_fusion.py --cff $cff --gene_bed $gene_bed --ref_fa $genome_fasta > $outdir/$(basename $cff).reann
  else 
    echo Annotate cff, no extraction of sequence surrounding breakpoint
    python reann_cff_fusion.py --cff $cff --gene_bed $gene_bed > $outdir/$(basename $cff).reann.NOSEQ
  fi
fi
# Assign .cff based on SEQ or NOSEQ
if [ $genome_fasta ]; then 
  cff=$outdir/$(basename $cff).reann
  echo cff $cff
else
  cff=$outdir/$(basename $cff).reann.NOSEQ
  echo cff $cff
fi

#Merge
cluster=$outdir/$(basename $cff).cluster
if [ $merge -eq 1 ]; then
  echo Merge cff
  bash RUN_cluster_genes_breakpoints.sh $cff $outdir $fusiontools > $cluster
fi

#output ANC_RT_SG file
if [ $output_ANC_RT_SG -eq 1 ]; then
  echo output ANC_RT_SG file
  python output_ANC_RT_SG.py $cluster > $cluster.ANC_RT_SG 
fi

#ReadThrough Callerfilter
if [ $RT_call_filter -eq 1 ]; then
  echo ReadThrough, callerfilter
  cat $cluster | grep ReadThrough > $outdir/$(basename $cluster).ReadThrough
  python callerfilter_num.py --cluster $cluster  --num_tools $num_tools | grep -v ReadThrough  > $outdir/$(basename $cluster).RT_filter.callerfilter.$num_tools
fi
cluster_RT_call=$outdir/$(basename $cluster).RT_filter.callerfilter.$num_tools 

# Blacklist Filter
if [ $blck_filter -eq 1 ]; then
  echo blacklist filter
  #blck_script_dir=/hpf/largeprojects/ccmbio/mapostolides/MODULES/FusionAnnotator/TEST_FusionAnnotator
#$blck_script_dir/blacklist_filter_recurrent_breakpoints.sh $cff $cluster_RT_call $outdir  > $outdir/$(basename $cluster).RT_filter.callerfilter.$num_tools.blck_filter
bash blacklist_filter_recurrent_breakpoints.sh $cff $cluster_RT_call $outdir $recurrent_bedpe > $outdir/$(basename $cluster).RT_filter.callerfilter.$num_tools.blck_filter
fi
cluster=$outdir/$(basename $cluster).RT_filter.callerfilter.$num_tools.blck_filter


# Adjacent Noncoding filter 
if [ $ANC_filter -eq 1 ]; then
  echo ANC adjacent noncoding filter
  python filter_adjacent_noncoding.py $cluster > $outdir/$(basename $cluster).ANC_filter  
fi
cluster=$outdir/$(basename $cluster).ANC_filter

#Rank and generate final.cluster
if [ $rank -eq 1 ]; then
   echo Rank and generate final.cluster 
  python rank_cluster_file.py $cluster > $outdir/final.cluster
fi
cluster=$outdir/final.cluster

#Benchmark
if [ $benchmark -eq 1 ]; then
   echo benchmark
  #/hpf/largeprojects/ccmbio/mapostolides/MODULES/RUN_BENCHMARKING_TOOLKIT/benchmarking_cluster-GENAP.sh $outdir $truth_set $cff $cluster true 
  benchmark_scripts=$fusiontools/FusionBenchmarking
  #fusionAnnotator=/hpf/tools/centos6/star-fusion/1.6.0/FusionAnnotator
  fusionAnnotator_dir=$fusiontools/FusionAnnotator
  if [ $FA -eq 1 ]; then
    #bash benchmarking_cluster-GENAP.sh $outdir $truth_set $cff $cluster $benchmark_scripts $fusionAnnotator_dir
    bash benchmarking_cluster-GENAP.sh $outdir $truth_set $cff $cluster $fusiontools FusionAnnotator
  else
    bash benchmarking_cluster-GENAP.sh $outdir $truth_set $cff $cluster $fusiontools 
  fi
fi


