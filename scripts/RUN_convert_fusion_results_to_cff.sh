#!/bin/bash
#convert_fusion_results_to_cff.py sample disease_name sample_type tool fusion_results_file out_dir

#raw_file_dir=/Users/mapostolides/MetaFusion/test_data/raw_output_STX16--RAE1
dataset=DREAM.SIM45.SIM52
raw_file_dir=/MetaFusion/test_data/caller_output_files/$dataset/star_seqr

outdir=/MetaFusion/RUNS/star_seqr.$dataset.cff_convert-Aug-19-2020
mkdir $outdir
disease=VALIDATION
sample_type=Tumor
tool=star_seqr
echo generating cff for $tool
samples=$(echo smc_rna_sim45 smc_rna_sim52)
for sample in ${samples[@]};do 
  result_file=$raw_file_dir/$sample\_STAR-SEQR/$sample\_STAR-SEQR_candidates.txt
  python convert_fusion_results_to_cff.py \
    $sample \
    $disease \
    $sample_type \
    $tool \
    $result_file \
    $outdir
done
