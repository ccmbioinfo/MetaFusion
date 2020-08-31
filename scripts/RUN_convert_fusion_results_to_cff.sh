#!/bin/bash
#convert_fusion_results_to_cff.py sample disease_name sample_type tool fusion_results_file out_dir

#raw_file_dir=/Users/mapostolides/MetaFusion/test_data/raw_output_STX16--RAE1
dataset=MELANOMA
raw_file_dir=/MetaFusion/test_data/caller_output_files/$dataset/star_seqr
date=Aug-21-2020

outdir=/MetaFusion/RUNS/star_seqr.$dataset.cff_convert-$date
mkdir $outdir
disease=melanoma_cml
sample_type=Tumor
tool=star_seqr
echo generating cff for $tool
samples=$(echo SRR018259 SRR018260 SRR018261 SRR018265 SRR018266 SRR018267 SRR018268 SRR018269)
for sample in ${samples[@]};do 
  #/MetaFusion/test_data/caller_output_files/MELANOMA/star_seqr/SRR018259/out_STAR-SEQR
  result_file=$raw_file_dir/$sample/out_STAR-SEQR/out_STAR-SEQR_candidates.txt
  python convert_fusion_results_to_cff.py \
    $sample \
    $disease \
    $sample_type \
    $tool \
    $result_file \
    $outdir
done
