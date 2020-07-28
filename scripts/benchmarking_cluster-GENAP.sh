#!/bin/bash
#INPUTS

#/hpf/largeprojects/ccmbio/mapostolides/MODULES/RUN_BENCHMARKING_TOOLKIT/benchmarking_cluster.sh $outdir $truth_fusions $cff $cluster $blacklist_filter $callerfilter2 $normal_filter  
outdir=$1
mkdir -p $outdir
truth_fusions=$2
cff=$3
cluster=$4
FUSION_BENCHMARK=$5
FUSION_ANNOTATOR=$6
#normal_filter=$5

#MODULE PATHS
#FUSION_ANNOTATOR=/hpf/tools/centos6/star-fusion/1.6.0/FusionAnnotator
perl_bin=/home/mapostolides/perl5/perlbrew/perls/perl-5.28.0/bin/perl
genome_lib_dir=/hpf/largeprojects/ccmbio/mapostolides/validate_fusion/test_star_star-fusion/GRCh37_v19_CTAT_lib_Feb092018.plug-n-play/ctat_genome_lib_build_dir

script_dir=/hpf/largeprojects/ccmbio/mapostolides/MODULES/FusionAnnotator/TEST_FusionAnnotator

#PUT CLUSTER INTO PROPER FORMAT
outfile=$outdir/cluster.preds.collected
echo -e "sample\tprog\tfusion\tJ\tS\tFID" > $outfile
#cat $cluster | awk '{FS=OFS="\t"}{print $15,"metacaller",$2"--"$3, $4, $5, $NF}' >> $outfile
# splits each call into one line per sample
while read -r line;do samples=$(echo $line | awk '{print $15}' | sed 's/,/ /g'); for sample in ${samples[@]};do echo $line | awk -v samp="$sample" '{$15=samp;print}' | sed 's/ /\t/g' | awk '{FS=OFS="\t"}{print $15,"metacaller",$2"--"$3, $4, $5, $NF}' ;done done < $cluster >> $outfile 

#/hpf/largeprojects/ccmbio/mapostolides/MODULES/FusionBenchmarking/resources/genes.aliases
echo Mapping gene partners to Gencode v19 genes 
   #${FUSION_BENCHMARK}/resources/genes.coords.gz \
#$perl_bin ${FUSION_BENCHMARK}/benchmarking/map_gene_symbols_to_gencode_FID.pl \
perl ${FUSION_BENCHMARK}/benchmarking/map_gene_symbols_to_gencode_FID.pl \
   $outfile \
   ${FUSION_BENCHMARK}/resources/genes.coords.HGNC_renamed_added.gz \
   ${FUSION_BENCHMARK}/resources/genes.aliases \
   > $outdir/$(basename $outfile).gencode_mapped
outfile=$outdir/$(basename $outfile).gencode_mapped


echo RUN FusionAnnotator
# "--full" parameter adds more detailed info to annotation
#$perl_bin ${FUSION_ANNOTATOR}/FusionAnnotator --annotate $outfile --genome_lib_dir $genome_lib_dir  -C 2 --full > $outdir/$(basename $outfile).wAnnot
perl ${FUSION_ANNOTATOR}/FusionAnnotator --annotate $outfile --genome_lib_dir $genome_lib_dir  -C 2 --full > $outdir/$(basename $outfile).wAnnot
outfile=$outdir/$(basename $outfile).wAnnot
exit 0
#SELECT FOR FUSIONS IN CANCER
#TRY ONLY THOSE CONFIRMED CANCER FUSIONS AND EXCLUDE "Individual genes of cancer relevance, which may show up in fusions" DATABASES
#Oncogene, ArcherDX_panel, FoundationOne_panel, OncocartaV1_panel, OncomapV4_panel
cat $outfile | grep 'FA_CancerSupp\|Mitelman\|chimerdb_omim\|chimerdb_pubmed\|ChimerKB\|ChimerPub\|ChimerSeq\|Cosmic\|YOSHIHARA_TCGA\|Klijn_CellLines\|Larsson_TCGA\|CCLE\|HaasMedCancer\|GUO2018CR_TCGA\|TumorFusionsNAR2018\|TCGA_StarF2019\|CCLE_StarF2019' > $outdir/$(basename $outfile).CANCER_FUSIONS
ids=$(cat $outdir/$(basename $outfile).CANCER_FUSIONS | cut -f 8| sort | uniq)
rm $outdir/$(basename $cluster).CANCER_FUSIONS
for id in ${ids[@]}; do
  cat $cluster | grep $id >> $outdir/$(basename $cluster).CANCER_FUSIONS;
done
cancer_cluster=$outdir/$(basename $cluster).CANCER_FUSIONS

#filter out anything which might be considered a "normal" tissue
#if [ $normal_filter == true ]; then echo normal_filter; 
#  cat $outfile | grep -v 'ConjoinG\|Babiceanu_Normal\|Greger_Normal\|HGNC_GENEFAM\|DGD_PARALOGS\|BodyMap\|GTEx' > $outdir/$(basename $outfile).NORMAL_FILTER 
#  cat $outfile | grep 'ConjoinG\|Babiceanu_Normal\|Greger_Normal\|HGNC_GENEFAM\|DGD_PARALOGS\|BodyMap\|GTEx' > $outdir/$(basename $outfile).NORMALS
#  outfile=$outdir/$(basename $outfile).NORMAL_FILTER
#fi 

# grep normal fusions, not a filter
cat $outfile | grep 'ConjoinG\|Babiceanu_Normal\|Greger_Normal\|HGNC_GENEFAM\|DGD_PARALOGS\|BodyMap\|GTEx' > $outdir/$(basename $outfile).NORMALS


echo Scoring of fusion predictions
${FUSION_BENCHMARK}/benchmarking/fusion_preds_to_TP_FP_FN_FID.pl --truth_fusions $truth_fusions --fusion_preds $outfile  --allow_reverse_fusion > $outdir/$(basename $outfile).scored 
outfile=$outdir/$(basename $outfile).scored

#GET TP/FP/FN
#echo -e TP\\tFP\\tFN > $outdir/$(basename $outfile).TP_FP_counts
#TP=$(cat $outfile | awk '$1=="TP"' | wc | awk '{print $1}')
#FP=$(cat $outfile | awk '$1=="FP"' | wc | awk '{print $1}')
#FN=$(cat $outfile | awk '$1=="FN"' | wc | awk '{print $1}')
#echo -e $TP\\t$FP\\t$FN >> $outdir/$(basename $outfile).TP_FP_counts

#TP
rm $outdir/$(basename $cluster).TP
ids=$( cat $outfile | awk '$1=="TP" || $1=="NA-TP"' | awk '{print $NF}' | sort | uniq)
touch $outdir/$(basename $cluster).TP 
for id in ${ids[@]}; do
  cat $cluster | grep $id >> $outdir/$(basename $cluster).TP ;
done
cat $outfile | awk '$1=="TP" || $1=="NA-TP"' > $outdir/$(basename $outfile).TP

#FP
rm $outdir/$(basename $cluster).FP
ids=$( cat $outfile | awk '$1=="FP" || $1=="NA-FP"' | awk '{print $NF}' | sort | uniq)
touch $outdir/$(basename $cluster).FP
for id in ${ids[@]}; do
  cat $cluster | grep $id >> $outdir/$(basename $cluster).FP ;
done
cat $outfile | awk '$1=="FP" || $1=="NA-FP"' > $outdir/$(basename $outfile).FP
#FN
cat $outfile | awk '$1=="FN"'  > $outdir/$(basename $outfile).FN

#GET TP/FP/FN 
# --> use cluster as FP to count merged calls and not per-sample
# --> use outfile as TP since this will allow for numerical match with truth set 
echo -e TP\\tFP\\tFN > $outdir/$(basename $cluster).TP_FP_counts
TP=$(cat $outfile | awk '$1=="TP"' | wc | awk '{print $1}')
FP=$(cat $outdir/$(basename $cluster).FP | wc | awk '{print $1}')
FN=$(cat $outfile | awk '$1=="FN"' | wc | awk '{print $1}')
#echo -e $TP\\t$FP\\t$FN >> $outdir/$(basename $outfile).TP_FP_counts
echo -e $TP\\t$FP\\t$FN >> $outdir/$(basename $cluster).TP_FP_counts
