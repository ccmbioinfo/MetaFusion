#!/bin/bash
caller_file_dir=/juno/work/ccs/noronhaa/tools/test_nf-core/anoronh4-forte_tests/full_test/results_review3/analysis

sampleinfo=/juno/work/ccs/pintoa1/fusion_report/metafusion/MetaFusion/testing_metafusion_sample_info.txt

outdir=/juno/work/ccs/pintoa1/fusion_report/metafusion/MetaFusion/testing_twelve_cff
mkdir $outdir
# Next, the tools used should be specified. Note that the names of the tools must exactly match "tool_name" at "/MetaFusion/caller_output_files/dataset/tool_name"
tools=$(echo arriba fusioncatcher starfusion)


while IFS= read -r sample_infor; do
  sample=`echo $sample_infor | awk '{print $1}'`
  disease=`echo $sample_infor | awk '{print $2}'`
  type=`echo $sample_infor | awk '{print $3}'`
  echo generating cff for $sample

	for tool in ${tools[@]};do
    echo $tool
    raw_file_dir=$caller_file_dir/$sample/$tool
    if [[ $tool = "arriba" ]]
    then
      result_file=$(ls $raw_file_dir/*.fusions.tsv)
    elif [[ $tool = "fusioncatcher" ]]
    then
      result_file=$(ls $raw_file_dir/*.fusion-genes.txt)
    else
      result_file=$(ls $raw_file_dir/*.fusion_predictions.tsv)
    fi
    echo  $sample, $tool, $result_file, $outdir, $disease, $type
    python convert_fusion_results_to_cff.py --sample $sample --disease_name $disease --sample_type $type --tool $tool --fusion_result_file $result_file --out_dir $outdir
  done

done < "$sampleinfo"

# Merge all .cff files into one combined file, "merged.cff"

files=$(ls $outdir/*cff | grep -v merged.cff)
cat $files > $outdir/merged.cff
