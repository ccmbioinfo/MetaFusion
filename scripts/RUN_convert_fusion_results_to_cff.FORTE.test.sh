#!/bin/bash

# This script converts a dataset, with result files from multiple callers, into a single "merged.cff" file.
# In order for this script to work correctly, the result files for each caller must be placed in a particular directory structure, as follows:

# /MetaFusion/path_to/caller_output_files/dataset/tool_name/sample_name/result_file.tsv
# A visual example of this directory structure can be found on the github wiki (https://github.com/ccmbioinfo/MetaFusion/wiki/How-to-generate-a-cff-file)

# The directory "/MetaFusion/caller_output_files" is the top-level directory containing result files, which is set here:
caller_file_dir=/juno/work/ccs/noronhaa/tools/test_nf-core/anoronh4-forte_tests/full_test/results_review3/analysis

# The sampleinfo file contains needed metadata to add to calls. The column headers must exactly match those specified below. An example of the TAB separated file is as follows:
#sample  disease sample_type
#BT474   BRCA      TP
#KPL4    BRCA      TP
#MCF7    BRCA      TP
#SKBR3   BRCA      TP
#In the above example, "TP" is used for Tumor tissue and "NT" is for normal tissue
#sampleinfo=/MetaFusion/Sandbox/BRCA_test/sampleinfo
sampleinfo=/juno/work/ccs/pintoa1/fusion_report/metafusion/MetaFusion/testing_metafusion_sample_info.txt

# The dataset variable is specified here. This string must exactly match the dataset name (i.e. the directory "dataset" in /MetaFusion/caller_output_files/dataset) as well as the "disease" field of the "sampleinfo file"
dataset=TheTwelve

# Next, the output directory, where the CFF files are generated, must be specified. It should be a directory different from where the caller output files are specified
outdir=/juno/work/ccs/pintoa1/fusion_report/metafusion/MetaFusion/testing_twelve_cff

# Next, the tools used should be specified. Note that the names of the tools must exactly match "tool_name" at "/MetaFusion/caller_output_files/dataset/tool_name"
tools=$(echo arriba fusioncatcher starfusion)

# The "sample_name" in "/MetaFusion/caller_output_files/dataset/tool_name/sample_name" will be used to name the sample, corresponds to the "sample" column in the sampleinfo file, and will appear in the output file as a field

# One last important note: There must be exactly and only ONE file inside each "/MetaFusion/caller_output_files/dataset/tool_name/sample_name/" directory. That single file is the fusion caller output file. Since often caller output file names differ depending on the user of the tool, this allows output files to be recognized regardless of how the file is named




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
    echo $sample, $tool, $result_file, $outdir, $disease, $type

    echo "python convert_fusion_results_to_cff.py \
      --sample $sample \
      --disease_name $disease \
      --sample_type $type \
      --tool $tool \
      --fusion_result_file $result_file \
      --outdir $outdir "
  done

done < "$sampleinfo"

# Merge all .cff files into one combined file, "merged.cff"

files=$(ls $outdir/*cff | grep -v merged.cff)
cat $files > $outdir/merged.cff
