#!/bin/bash
#outdir=$1
#mkdir -p $outdir
#cff=$2
#gene_bed=$3
#truth_fusions=$4
#num_tools=$5
#gene_info=$6
#genome_fasta=$7

#ARGS
OTHER_ARGUMENTS=()

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
        *)
        OTHER_ARGUMENTS+=("$1")
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
#
rename=1
annotate=1
merge=1
output_ANC_RT_SG=1
RT_call_filter=1
blck_filter=1
ANC_filter=1
benchmark=1

fusiontools=/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin
#Check CFF file format:
#Remove entries with nonconformming chromosome name
cat $cff | awk '$1 ~ /[0-9XY]/ && $4 ~ /[0-9XY]/ ' > $outdir/$(basename $cff).reformat 
cff=$outdir/$(basename $cff).reformat
#NEED TO INSERT +/-/NA for strand, make NA if other

#Rename cff
if [ $rename -eq 1 ]; then
  echo Rename cff
  python $fusiontools/rename_cff_file_genes-GENAP.py $cff $gene_info > $outdir/$(basename $cff).renamed
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
  /hpf/largeprojects/ccmbio/mapostolides/MODULES/RUN_BENCHMARKING_TOOLKIT/benchmarking_cluster-GENAP.sh $outdir $truth_set $cff $cluster true 
fi


