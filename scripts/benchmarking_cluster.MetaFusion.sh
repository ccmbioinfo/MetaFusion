#!/bin/bash
outdir=$1
mkdir -p $outdir
truth_fusions=$2
cluster=$3
fusiontools=$4
FusionAnnotator=$5

FUSION_BENCHMARK=$fusiontools/FusionBenchmarking
FUSION_ANNOTATOR=$fusiontools/FusionAnnotator
#normal_filter=$5

#RUN FusionAnnotator if path to "FUSION_ANNOTATOR" is specified
if [ $FusionAnnotator ]; then
    FA=1
else
	FA=0
fi

#MODULE PATHS
genome_lib_dir=$(dirname $fusiontools)/reference_files/ctat_genome_lib_build_dir

#PUT CLUSTER INTO PROPER FORMAT
outfile=$outdir/cluster.preds.collected
echo -e "sample\tprog\tfusion\tJ\tS\tFID" > $outfile
# splits each call into one line per sample
# CLUSTER FORMAT
#while read -r line;do samples=$(echo $line | awk '{print $15}' | sed 's/,/ /g'); for sample in ${samples[@]};do echo $line | awk -v samp="$sample" '{$15=samp;print}' | sed 's/ /\t/g' | awk '{FS=OFS="\t"}{print $15,"metacaller",$2"--"$3, $4, $5, $NF}' ;done done < $cluster  | grep -v fusion_IDs  >> $outfile 
#CLUSTER SUBSET FORMAT
while read -r line;do samples=$(echo $line | awk '{print $13}' | sed 's/,/ /g'); for sample in ${samples[@]};do echo $line | awk -v samp="$sample" '{$13=samp;print}' | sed 's/ /\t/g' | awk '{FS=OFS="\t"}{print $13, "metacaller", $1"--"$2, $7, $8, $NF}' ;done done < $cluster  | grep -v fusion_IDs  >> $outfile 

#$13, "metacaller", ,$1"--"$2, $7, $8, $16, 
#SRR018266   metacaller  AAMDC,RP11-91P24.6--SLC12A7 62  56  F00000029,F00000942,F00000690,F00001130,F00000047,F00000499,F00000285,F00000284,F00000554,F00000498,F00000046

echo Mapping gene partners to Gencode v19 genes 
perl ${FUSION_BENCHMARK}/benchmarking/map_gene_symbols_to_gencode_FID.pl \
   $outfile \
   ${FUSION_BENCHMARK}/resources/genes.coords.HGNC_renamed_added.gz \
   ${FUSION_BENCHMARK}/resources/genes.aliases \
   > $outdir/$(basename $outfile).gencode_mapped
outfile=$outdir/$(basename $outfile).gencode_mapped

#BEGIN FusionAnnotator
#if [ $FA -eq 1 ]; then

#echo RUN FusionAnnotator
## "--full" parameter adds more detailed info to annotation
#perl ${FUSION_ANNOTATOR}/FusionAnnotator --annotate $outfile --genome_lib_dir $genome_lib_dir  -C 2 --full > $outdir/$(basename $outfile).wAnnot
#outfile=$outdir/$(basename $outfile).wAnnot
#
##SELECT FOR FUSIONS IN CANCER, ONLY THOSE CONFIRMED CANCER FUSIONS AND EXCLUDE "Individual genes of cancer relevance, which may show up in fusions" DATABASES
##Oncogene, ArcherDX_panel, FoundationOne_panel, OncocartaV1_panel, OncomapV4_panel
#cat $outfile | grep 'FA_CancerSupp\|Mitelman\|chimerdb_omim\|chimerdb_pubmed\|ChimerKB\|ChimerPub\|ChimerSeq\|Cosmic\|YOSHIHARA_TCGA\|Klijn_CellLines\|Larsson_TCGA\|CCLE\|HaasMedCancer\|GUO2018CR_TCGA\|TumorFusionsNAR2018\|TCGA_StarF2019\|CCLE_StarF2019' > $outdir/$(basename $outfile).CANCER_FUSIONS
#ids=$(cat $outdir/$(basename $outfile).CANCER_FUSIONS | cut -f 8| sort | uniq)
#
##add header
#echo \#cluster_type gene1 gene2 max_split_cnt max_span_cnt sample_type disease tools inferred_fusion_type gene1_on_bnd gene1_close_to_bnd gene2_on_bnd gene2_close_to_bnd dna_supp samples chr1 breakpoint_1 chr2 breakpoint_2 cancer_db_hits captured_reads_normal_mean fusion_IDs | sed 's/ /\t/g' > $outdir/$(basename $cluster).CANCER_FUSIONS 
#
#for id in ${ids[@]}; do
#  cat $cluster | grep $id >> $outdir/$(basename $cluster).CANCER_FUSIONS;
#done
#cancer_cluster=$outdir/$(basename $cluster).CANCER_FUSIONS

# NORMAL FILTER
#filter out anything which might be considered a "normal" tissue
#if [ $normal_filter == true ]; then echo normal_filter; 
#  cat $outfile | grep -v 'ConjoinG\|Babiceanu_Normal\|Greger_Normal\|HGNC_GENEFAM\|DGD_PARALOGS\|BodyMap\|GTEx' > $outdir/$(basename $outfile).NORMAL_FILTER 
#  cat $outfile | grep 'ConjoinG\|Babiceanu_Normal\|Greger_Normal\|HGNC_GENEFAM\|DGD_PARALOGS\|BodyMap\|GTEx' > $outdir/$(basename $outfile).NORMALS
#  outfile=$outdir/$(basename $outfile).NORMAL_FILTER
#fi 

# grep normal fusions, not a filter
#cat $outfile | grep 'ConjoinG\|Babiceanu_Normal\|Greger_Normal\|HGNC_GENEFAM\|DGD_PARALOGS\|BodyMap\|GTEx' > $outdir/$(basename $outfile).NORMALS

#END FusionAnnotator
#fi

echo Scoring of fusion predictions
${FUSION_BENCHMARK}/benchmarking/fusion_preds_to_TP_FP_FN_FID.pl --truth_fusions $truth_fusions --fusion_preds $outfile  --allow_reverse_fusion > $outdir/$(basename $outfile).scored 
outfile=$outdir/$(basename $outfile).scored

#TP
#generate TP cluster file
#add header
echo \#gene1  gene2   chr1    breakpoint_1    chr2    breakpoint_2    max_split_cnt   max_span_cnt    sample_type disease tools   inferred_fusion_type    samples cancer_db_hits  fusion_IDs | sed 's/ \+/\t/g' > $outdir/$(basename $cluster).TP
ids=$( cat $outfile | awk '$1=="TP" || $1=="NA-TP"' | awk '{print $NF}' | sort | uniq)
for id in ${ids[@]}; do
  cat $cluster | grep $id >> $outdir/$(basename $cluster).TP ;
done

#generate TP outfile in benchmarking toolkit format
cat $outfile | awk '$1=="TP" || $1=="NA-TP"' > $outdir/$(basename $outfile).TP

#FP
ids=$( cat $outfile | awk '$1=="FP" || $1=="NA-FP"' | awk '{print $NF}' | sort | uniq)
#generate FP cluster file
#add header
#echo \#cluster_type gene1 gene2 max_split_cnt max_span_cnt sample_type disease tools inferred_fusion_type gene1_on_bnd gene1_close_to_bnd gene2_on_bnd gene2_close_to_bnd dna_supp samples chr1 breakpoint_1 chr2 breakpoint_2 cancer_db_hits captured_reads_normal_mean fusion_IDs | sed 's/ /\t/g' > $outdir/$(basename $cluster).FP
echo \#gene1  gene2   chr1    breakpoint_1    chr2    breakpoint_2    max_split_cnt   max_span_cnt    sample_type disease tools   inferred_fusion_type    samples cancer_db_hits  fusion_IDs | sed 's/ \+/\t/g' > $outdir/$(basename $cluster).FP
for id in ${ids[@]}; do
  cat $cluster | grep $id >> $outdir/$(basename $cluster).FP ;
done

#generate FP outfile in benchmarking toolkit format
cat $outfile | awk '$1=="FP" || $1=="NA-FP"' > $outdir/$(basename $outfile).FP

#FN
#generate FN outfile in benchmarking toolkit format
cat $outfile | awk '$1=="FN"'  > $outdir/$(basename $outfile).FN

#GET TP/FP/FN counts 
# --> use outfile as TP since this will allow for numerical match with truth set 
echo -e TP\\tFP\\tFN > $outdir/$(basename $cluster).TP_FP_counts
TP=$(cat $outfile | awk '$1=="TP"' | wc | awk '{print $1}')
FP=$(cat $outdir/$(basename $outfile).FP | wc | awk '{print $1}')
FN=$(cat $outfile | awk '$1=="FN"' | wc | awk '{print $1}')
echo -e $TP\\t$FP\\t$FN >> $outdir/$(basename $cluster).TP_FP_counts
