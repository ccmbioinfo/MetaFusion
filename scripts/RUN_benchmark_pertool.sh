#!/bin/bash

#Change date to current date
date=Sept-2-2020

#DATASETS
sim45_sim52=0
brca_4=0
beers_neg=0
sim_50=0
sim101=0
melanoma=1


#REFERENCE FILES 
gene_bed=/MetaFusion/reference_files/ens_known_genes.renamed.ENSG.bed
gene_info=/MetaFusion/reference_files/Homo_sapiens.gene_info
genome_fasta=/MetaFusion/reference_files/human_g1k_v37_decoy.fasta
recurrent_bedpe=/MetaFusion/reference_files/blacklist_breakpoints.bedpe
fusiontools=/MetaFusion/scripts
runs_dir=/MetaFusion/RUNS


# SIM45.SIM52.combined
if [ $sim45_sim52 -eq 1 ]; then
echo SIM45.SIM52
outdir=$runs_dir/SIM45.SIM52.pertool_benchmark.$date
echo generating output in $outdir
mkdir $outdir
cff=/MetaFusion/test_data/cff/dream.sim45.sim52.cff
truth_fusions=/MetaFusion/test_data/truth_sets/dream.sim45.sim52.truth_set.dat

bash $fusiontools/benchmarking_cff_pertool.sh $outdir $truth_fusions $cff $fusiontools

fi

#BT474.KPL4.MCF7.SKBR3
if [ $brca_4 -eq 1 ]; then
echo BT474.KPL4.MCF7.SKBR3
outdir=$runs_dir/BT474.KPL4.MCF7.SKBR3.pertool_benchmark.$date
echo generating output in $outdir
truth_fusions=/MetaFusion/test_data/truth_sets/BRCA.truth_set.dat
cff=/MetaFusion/test_data/cff/BRCA.cff

bash $fusiontools/benchmarking_cff_pertool.sh $outdir $truth_fusions $cff $fusiontools
fi

#NEGATIVE CONTROL BEERS
if [ $beers_neg -eq 1 ]; then
echo NEGATIVE CONTROL BEERS
outdir=$runs_dir/NEG_CONTROL_BEERS.pertool_benchmark.$date
echo generating output in $outdir
cff=/MetaFusion/test_data/cff/beers_neg.cff 
truth_fusions=/MetaFusion/test_data/truth_sets/BRCA.truth_set.dat

bash $fusiontools/benchmarking_cff_pertool.sh $outdir $truth_fusions $cff $fusiontools
fi

# SIM50 2500 fusions files:
if [ $sim_50 -eq 1 ]; then
echo SIM50
outdir=$runs_dir/SIM50.pertool_benchmark.$date
echo generating output in $outdir
cff=/MetaFusion/test_data/cff/sim50.cff
truth_fusions=/MetaFusion/test_data/truth_sets/sim50.truth_set.dat

bash $fusiontools/benchmarking_cff_pertool.sh $outdir $truth_fusions $cff $fusiontools
fi

#SIM101 2500 fusions files, same truth set as SIM50
if [ $sim101 -eq 1 ]; then
echo SIM101
outdir=$runs_dir/SIM101.pertool_benchmark.$date
echo generating output in $outdir
cff=/MetaFusion/test_data/cff/sim101.cff
truth_fusions=/MetaFusion/test_data/truth_sets/sim101.truth_set.dat

bash $fusiontools/benchmarking_cff_pertool.sh $outdir $truth_fusions $cff $fusiontools
fi

## Melanoma and CML
if [ $melanoma -eq 1 ]; then
echo MELANOMA and CML
outdir=$runs_dir/melanoma.CML.pertool_benchmark.$date
#outdir=$runs_dir/melanoma.CML.pertool_benchmark.truth_set.NO_DUPS.$date
echo generating output in $outdir
cff=/MetaFusion/test_data/cff/melanoma.cff
#cff=/MetaFusion/test_data/cff/melanoma.no_SRR018269.cff
truth_fusions=/MetaFusion/test_data/truth_sets/melanoma.truth_set.dat
#truth_fusions=/MetaFusion/test_data/truth_sets/melanoma.truth_set.no_SRR018269.dat

bash $fusiontools/benchmarking_cff_pertool.sh $outdir $truth_fusions $cff $fusiontools
fi
