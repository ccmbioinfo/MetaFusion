#!/bin/bash

cff=$1
cluster=$2
outdir=$3
recurrent_bedpe=$4

#BLOCKLIST FILTER
cff_bedpe=$(basename $cff).bedpe
cat $cff | awk '{FS = OFS = "\t"}{print $1,$2,$2+1,$4,$5,$5+1,$31}' | sed 's/chr//g' > $outdir/$cff_bedpe 
pairToPair -a $recurrent_bedpe -b $outdir/$cff_bedpe -slop 5 > $outdir/recurrent_filtered_fusions.bedpe

# Make sorted blocklist file
ids=$(cat $outdir/recurrent_filtered_fusions.bedpe | awk '{print $NF}' | sort | uniq)
cat $outdir/recurrent_filtered_fusions.bedpe | awk '{print $NF}' | sort | uniq > $outdir/FIDS
for FID in ${ids[@]};do cat $cluster | grep $FID ;done | sort | uniq  >  $outdir/$(basename $cluster).BLOCKLIST

# Sort cluster, required for "comm" command
cat $cluster | sort > $outdir/$(basename $cluster).sorted

# Remove lines of BLOCKLIST file from cluster, creating .blck_filter file
comm -2 -3 $outdir/$(basename $cluster).sorted $outdir/$(basename $cluster).BLOCKLIST 
