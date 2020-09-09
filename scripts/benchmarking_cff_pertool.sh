#!/bin/bash
#INPUTS
outdir=$1
mkdir -p $outdir
truth_fusions=$2
cff=$3
fusiontools=$4
FUSION_BENCHMARK=$fusiontools/FusionBenchmarking

#CALLERS
#callers=$(echo arriba star_fusion star_seqr defuse ericscript integrate fusionmap chimerascan)
callers=$(echo arriba star_fusion star_seqr defuse ericscript integrate fusionmap)

# Add FID to unannotated .cff
python $fusiontools/add_FID_to_cff.py $cff > $outdir/$(basename $cff).FID
cff=$outdir/$(basename $cff).FID

#PUT CFF INTO PROPER FORMAT
outfile=$outdir/cff.preds.collected
echo -e "sample\tprog\tfusion\tJ\tS\tFID" > $outfile
#cat $cff | awk '{FS=OFS="\t"}{print $8,$11,$14"--"$16, $12, $13, $31}' >> $outfile
#cat $cff | sed 's/(.\+)//g' | sed 's/\//,/g' | awk '{FS=OFS="\t"}{print $8,$11,$14"--"$16, $12, $13, $31}' >> $outfile

echo Mapping gene partners to Gencode v19 genes 
$perl_bin ${FUSION_BENCHMARK}/benchmarking/map_gene_symbols_to_gencode_FID.pl \
   $outfile \
   ${FUSION_BENCHMARK}/resources/genes.coords.HGNC_renamed_added.gz \
   ${FUSION_BENCHMARK}/resources/genes.aliases \
   > $outdir/$(basename $outfile).gencode_mapped
outfile=$outdir/$(basename $outfile).gencode_mapped

echo Scoring of fusion predictions
${FUSION_BENCHMARK}/benchmarking/fusion_preds_to_TP_FP_FN_FID.pl --truth_fusions $truth_fusions --fusion_preds $outfile  --allow_reverse_fusion > $outdir/$(basename $outfile).scored 
outfile=$outdir/$(basename $outfile).scored


#GET TP/FP/FN
echo -e caller\\tTP\\tFP\\tFN > $outdir/$(basename $outfile).TP_FP_counts
for caller in ${callers[@]};do 
  # per-caller cff and outfile
  cat $cff | grep $caller > $cff.$caller 
  outfile_caller=$outdir/$(basename $outfile).$caller
  cat $outfile | grep $caller > $outfile_caller

  # TP/FP/FN outfiles
  cat $outfile_caller | awk '$1=="TP" || $1=="NA-TP"' > $outfile_caller.TP
  cat $outfile_caller | awk '$1=="FP" || $1=="NA-FP"' > $outfile_caller.FP
  cat $outfile_caller | awk '$1=="FN"' > $outfile_caller.FN
  # TP/FP/FN counts
  TP=$(cat $outfile_caller | awk '$1=="TP"'  | wc | awk '{print $1}') 
  FP=$(cat $outfile_caller.FP  | awk '$1=="FP"' | wc | awk '{print $1}')
  FN=$(cat $outfile_caller.FN | awk '$1=="FN"'| wc | awk '{print $1}')
  echo -e $caller\\t$TP\\t$FP\\t$FN >> $outdir/$(basename $outfile).TP_FP_counts

  #generate per-caller .cff for TP and FP
  #TP_IDs=$( cat $outfile_caller.TP | awk '{print $NF}' | sort | uniq)
  #FP_IDs=$( cat $outfile_caller.FP | awk '{print $NF}' | sort | uniq)
  ##TP
  #rm $cff.$caller.TP
  #for id in ${TP_IDs[@]}; do
  #  cat $cff.$caller | grep $id >> $cff.$caller.TP;
  #done
  ##FP
  #rm $cff.$caller.FP
  #for id in ${FP_IDs[@]}; do
  #  cat $cff.$caller | grep $id >> $cff.$caller.FP;
  #done
done

