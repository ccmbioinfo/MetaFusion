#!/bin/bash
#convert_fusion_results_to_cff.py sample disease_name sample_type tool fusion_results_file out_dir

#raw_file_dir=/Users/mapostolides/MetaFusion/test_data/raw_output_STX16--RAE1
raw_file_dir=/MetaFusion/test_data/raw_output_STX16--RAE1
outdir=/MetaFusion/RUNS/tmp
mkdir $outdir
callers=$(echo star_fusion fusionmap arriba integrate defuse ericscript star_seqr)
sample=BT474
disease=BRCA
sample_type=Tumor
echo ${callers[@]}
for caller in ${callers[@]};do 
  tool=$caller
  echo generating cff for $tool
  result_file=$raw_file_dir/$tool\_raw.STX16--RAE1
  python convert_fusion_results_to_cff.py \
    $sample \
    $disease \
    $sample_type \
    $tool \
    $result_file \
    $outdir
done 
exit 0
convert_fusion_results_to_cff.py \
  $sample \
  ../sampleinfo.BT474.KPL4.MCF7.SKBR3 \
  fusionmap \
  fusions/fusionmap/BT474/02_RNA.FusionReport.txt \
  fusions/cff && \
convert_fusion_results_to_cff.py \
  $sample \
  ../sampleinfo.BT474.KPL4.MCF7.SKBR3 \
  arriba \
  fusions/arriba/BT474/fusions.tsv \
  fusions/cff && \
convert_fusion_results_to_cff.py \
  $sample \
  ../sampleinfo.BT474.KPL4.MCF7.SKBR3 \
  integrate \
  fusions/integrate/BT474/breakpoints.cov.tsv \
  fusions/cff && \
convert_fusion_results_to_cff.py \
  $sample \
  ../sampleinfo.BT474.KPL4.MCF7.SKBR3 \
  star_seqr \
  fusions/star_seqr/BT474/out_STAR-SEQR/out_STAR-SEQR_candidates.txt \
  fusions/cff && \
convert_fusion_results_to_cff.py \
  $sample \
  ../sampleinfo.BT474.KPL4.MCF7.SKBR3 \
  defuse \
  fusions/defuse/BT474/results.filtered.tsv \
  fusions/cff && \
convert_fusion_results_to_cff.py \
  $sample \
  ../sampleinfo.BT474.KPL4.MCF7.SKBR3 \
  ericscript \
  fusions/ericscript/BT474/fusion.results.filtered.tsv \
  fusions/cff
