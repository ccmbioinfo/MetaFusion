#!/bin/bash

# This script converts a dataset, with result files from multiple callers, into a single "merged.cff" file.
# In order for this script to work correctly, the result files for each caller must be placed in a particular directory structure, as follows:

# /MetaFusion/path_to/caller_output_files/dataset/tool_name/sample_name/result_file.tsv
# A visual example of this directory structure can be found on the github wiki (https://github.com/ccmbioinfo/MetaFusion/wiki/How-to-generate-a-cff-file)

# The directory "/MetaFusion/caller_output_files" is the top-level directory containing result files, which is set here: 
caller_file_dir=/MetaFusion/test_data/caller_output_files/

# The sampleinfo file contains needed metadata to add to calls. The column headers must exactly match those specified below. An example of the TAB separated file is as follows:
#sample  disease sample_type
#BT474   BRCA      TP
#KPL4    BRCA      TP
#MCF7    BRCA      TP
#SKBR3   BRCA      TP
#In the above example, "TP" is used for Tumor tissue and "NT" is for normal tissue
#sampleinfo=/MetaFusion/Sandbox/BRCA_test/sampleinfo
sampleinfo=/MetaFusion/test_data/sampleinfo/sampleinfo.BRCA

# The dataset variable is specified here. This string must exactly match the dataset name (i.e. the directory "dataset" in /MetaFusion/caller_output_files/dataset) as well as the "disease" field of the "sampleinfo file" 
dataset=BRCA

# Next, the output directory, where the CFF files are generated, must be specified. It should be a directory different from where the caller output files are specified 
outdir=/MetaFusion/Sandbox/BRCA_cff

# Next, the tools used should be specified. Note that the names of the tools must exactly match "tool_name" at "/MetaFusion/caller_output_files/dataset/tool_name"
tools=$(echo arriba star_fusion star_seqr defuse ericscript integrate fusionmap)

# The "sample_name" in "/MetaFusion/caller_output_files/dataset/tool_name/sample_name" will be used to name the sample, corresponds to the "sample" column in the sampleinfo file, and will appear in the output file as a field

# One last important note: There must be exactly and only ONE file inside each "/MetaFusion/caller_output_files/dataset/tool_name/sample_name/" directory. That single file is the fusion caller output file. Since often caller output file names differ depending on the user of the tool, this allows output files to be recognized regardless of how the file is named

for tool in ${tools[@]};do 

raw_file_dir=$caller_file_dir/$dataset/$tool
result_files=$(ls $raw_file_dir/*/*)

echo generating cff for $tool
for result_file in ${result_files[@]};do
	sample=$(basename $(dirname $result_file))
	echo $sample, $tool, $result_file, $outdir
  python convert_fusion_results_to_cff.py \
    --sample $sample \
    --sample_info_file $sampleinfo \
    --tool $tool \
    --fusion_result_file $result_file \
    --outdir $outdir
done

done

# Merge all .cff files into one combined file, "merged.cff"

files=$(ls $outdir/*cff | grep -v merged.cff)
cat $files > $outdir/merged.cff
