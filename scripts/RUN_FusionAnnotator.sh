#!/bin/bash
outdir=$1
mkdir -p $outdir
cluster=$2
fusiontools=$3

FUSION_ANNOTATOR=$fusiontools/FusionAnnotator
FUSION_BENCHMARK=$fusiontools/FusionBenchmarking
#normal_filter=$5

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

#echo Mapping gene partners to Gencode v19 genes 
perl ${FUSION_BENCHMARK}/benchmarking/map_gene_symbols_to_gencode_FID.pl \
   $outfile \
   ${FUSION_BENCHMARK}/resources/genes.coords.HGNC_renamed_added.gz \
   ${FUSION_BENCHMARK}/resources/genes.aliases \
   > $outdir/$(basename $outfile).gencode_mapped
outfile=$outdir/$(basename $outfile).gencode_mapped

#BEGIN FusionAnnotator

#echo RUN FusionAnnotator
# "--full" parameter adds more detailed info to annotation
perl ${FUSION_ANNOTATOR}/FusionAnnotator --annotate $outfile --genome_lib_dir $genome_lib_dir  -C 2 --full > $outdir/$(basename $outfile).wAnnot
outfile=$outdir/$(basename $outfile).wAnnot

#SELECT FOR FUSIONS IN CANCER, ONLY THOSE CONFIRMED CANCER FUSIONS AND EXCLUDE "Individual genes of cancer relevance, which may show up in fusions" DATABASES
#Oncogene, ArcherDX_panel, FoundationOne_panel, OncocartaV1_panel, OncomapV4_panel
cat $outfile | grep 'FA_CancerSupp\|Mitelman\|chimerdb_omim\|chimerdb_pubmed\|ChimerKB\|ChimerPub\|ChimerSeq\|Cosmic\|YOSHIHARA_TCGA\|Klijn_CellLines\|Larsson_TCGA\|CCLE\|HaasMedCancer\|GUO2018CR_TCGA\|TumorFusionsNAR2018\|TCGA_StarF2019\|CCLE_StarF2019' > $outdir/$(basename $outfile).CANCER_FUSIONS

# grep normal fusions, not a filter
cat $outfile | grep 'ConjoinG\|Babiceanu_Normal\|Greger_Normal\|HGNC_GENEFAM\|DGD_PARALOGS\|BodyMap\|GTEx' > $outdir/$(basename $outfile).NORMALS
