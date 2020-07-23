#!/bin/bash

#module load bedtools

cff=$1
cluster=$2
outdir=$3

#BLACKLIST FILTER
recurrent_bedpe=/hpf/largeprojects/ccmbio/mapostolides/validate_fusion/stjude_validation/Analysis_manuscript/OUTDIR_APR-3-2020-sim45-52/recurrent_breakpoints_filtering/blacklist_breakpoints.bedpe
cff_bedpe=$(basename $cff).bedpe
cat $cff | awk '{FS = OFS = "\t"}{print $1,$2,$2+1,$4,$5,$5+1,$31}' | sed 's/chr//g' > $outdir/$cff_bedpe 
pairToPair -a $recurrent_bedpe -b $outdir/$cff_bedpe -slop 5 > $outdir/recurrent_filtered_fusions.bedpe

# Make sorted blacklist file
ids=$(cat $outdir/recurrent_filtered_fusions.bedpe | awk '{print $NF}' | sort | uniq)
for FID in ${ids[@]};do cat $cluster | grep $FID ;done | sort | uniq  >  $outdir/$(basename $cluster).BLACKLIST

# Sort cluster, required for "comm" command
cat $cluster | sort > $outdir/$(basename $cluster).sorted

# Remove lines of BLACKLIST file from cluster, creating .blck_filter file
comm -2 -3 $outdir/$(basename $cluster).sorted $outdir/$(basename $cluster).BLACKLIST #> $outdir/$(basename $cluster).blck_filter



#SCRAP
#cluster_norm_rm=$outdir/$(basename $cluster).NORMALS_REMOVED
#for FID in ${ids[@]};do cat $cluster | grep $FID ;done #| sort | uniq #>$outdir/$(basename $cluster).BLACKLIST
#cat $cluster | grep $grep_regex > $outdir/$(basename $cluster).BLACKLIST

#exit 0
# select for Tumor and 2+ callers
#cat $outdir/$(basename $cluster_norm_rm).BLACKLIST_FILTERED | awk '$6=="Tumor"' | awk '$8 ~ /.,./ {print $0}' > $outdir/$(basename $cluster_norm_rm).BLACKLIST_FILTERED.Tumor.callerfilter2
#
##select for more than one sample
#cat $outdir/$(basename $cluster_norm_rm).BLACKLIST_FILTERED | awk '$15 ~ /.,./ {print $0}' > $outdir/$(basename $cluster_norm_rm).BLACKLIST_FILTERED.2_plus_samples 
#

#grep_regex=$(echo $ids | sed 's/ /\\|/g' | sed 's/^/\x27/g' | sed 's/$/\x27/g')
#grep_regex=$(echo $ids | sed 's/ /|/g' | sed 's/^/\x27/g' | sed 's/$/\x27/g')
#echo $grep_regex
#cat $cluster_norm_rm | grep -v $grep_regex > $outdir/$(basename $cluster_norm_rm).BLACKLIST_FILTERED

#OUTPUT FILE TO SDTOUT
#echo $grep_regex
#cat $cluster | grep -Ev $grep_regex 
