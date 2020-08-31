#!/bin/bash
#NOTES
#$reann_test_dir/cff_files/NOTES-merged.cff_RPRD2--LAMC2
source setup.sh

date=JULY-27-2020
#CONSTANTS
reann_test_dir=/hpf/largeprojects/ccmbio/mapostolides/MODULES/MetaFusion/reann_cff_fusion_testing
gene_bed_total=/hpf/largeprojects/ccmbio/mapostolides/gene_fusion/annotation/ens_known_genes.renamed.ENSG.bed
bed_dir=$reann_test_dir/bed_files
cff_dir=$reann_test_dir/cff_files

#/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/reann_cff_fusion.py 
#/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/generate_common_fusion_stats.py 
run_reann_cluster (){
  outdir=$1
  mkdir -p $outdir
  cff=$2
  gene_bed=$3
  fusiontools=/hpf/largeprojects/ccmbio/mapostolides/MODULES/MetaFusion/scripts
  genome_fasta=/hpf/largeprojects/ccmbio/mapostolides/gene_fusion/pipeline/config_reference_files/human_g1k_v37_decoy.fasta
  #cat $cff
  python $fusiontools/reann_cff_fusion_exon_test2.py $cff $gene_bed $genome_fasta > $outdir/$(basename $cff).exons

  #cat $outdir/$(basename $cff).reann 
  #python $fusiontools/generate_common_fusion_stats.py  $outdir/$(basename $cff).reann > $outdir/$(basename $cff).reann.cluster 
}
#NOTCH1--NUP214
cff=$reann_test_dir/cff_files/NOTCH1--NUP214.cff
cff=$reann_test_dir/cff_files/NOTCH1--NUP214.134062675.cff
outdir=$reann_test_dir/outdir-$date
gene_bed=$reann_test_dir/ens_known_genes.NOTCH1--NUP214.bed
#gene_bed=$gene_bed_total
#run_reann_cluster $outdir $cff $gene_bed

#exit 0

#4BRCA
cff=/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/testing_pipeline/BT474.KPL4.MCF7.SKBR3.benchmark.July-29-2020/merged.cff.reformat.renamed.reann
outdir=$reann_test_dir/outdir-$date
gene_bed=$gene_bed_total
run_reann_cluster $outdir $cff $gene_bed

#STX16--RAE1
cff=/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/testing_pipeline/BT474.KPL4.MCF7.SKBR3.benchmark.June-18-2020/STX16-NPEPL1,STX16--RAE1,MIR5095.cff.renamed.reann-SAVE
outdir=$reann_test_dir/outdir-$date
gene_bed=$reann_test_dir/ens_known_genes.STX16--RAE1.renamed.bed
run_reann_cluster $outdir $cff $gene_bed

exit 0 
#outdir=$reann_test_dir/outdir.ERBIN--RASD.May-21-2020.gene_bed_subset 
outdir=$reann_test_dir/outdir-$date
#cff=$reann_test_dir/cff_files/ERBIN--RASD.cff.renamed
#cff=$reann_test_dir/cff_files/ERBIN--RASD.cff.renamed.2entries
cff=$reann_test_dir/cff_files/ERBIN--RASD.cff.renamed.1entry 
gene_bed=$reann_test_dir/ens_known_genes.ERBIN--RASD.renamed.bed
run_reann_cluster $outdir $cff $gene_bed
exit 0

# EFTUD1 reannotated to ELF1 incorrectly.
# /hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/reann_cff_fusion_testing/cff_files/TLCD2--EFTUD1.cff

#GeneFusion misannotated as TruncatedCoding -- Actually PVT1 is lncRNA, so annotation is correct!
outdir=$reann_test_dir/outdirs/outdir.PVT1--PDGFB.June1-2020.gene_bed_subset 
cff=/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/reann_cff_fusion_testing/cff_files/PVT1--PDGFB.cff.arriba.1entry
#gene_bed=$reann_test_dir/ens_known_genes.renamed.SLC5A1--PNMA6A.bed
gene_bed=$reann_test_dir/ens_known_genes.renamed.PVT1--PDGFB.ENSG.bed
run_reann_cluster $outdir $cff $gene_bed
exit 0

#GeneFusion misannotated as TruncatedCoding
outdir=$reann_test_dir/outdirs/outdir.SLC5A1--PNMA6A.May-25-2020.gene_bed_subset 
cff=/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/reann_cff_fusion_testing/cff_files/SLC5A1--PNMA6A_Trunc_miscat.cff
#gene_bed=$reann_test_dir/ens_known_genes.renamed.SLC5A1--PNMA6A.bed
gene_bed=$reann_test_dir/ens_known_genes.renamed.SLC5A1--PNMA6A.bed.1entry_removed
run_reann_cluster $outdir $cff $gene_bed
exit 0


#outdir=$reann_test_dir/outdir.DCPS--TIRAP.May-21-2020.gene_bed_subset
outdir=$reann_test_dir/outdir.DCPS--TIRAP_RP11-712L6.5.May-21-2020.gene_bed_subset
#DCPS--TIRAP- SameGene bug
cff=$reann_test_dir/DCPS--TIRAP/DCPS--TIRAP.cff
#gene_bed=$reann_test_dir/ens_known_genes_DCPS--TIRAP.bed
gene_bed=$reann_test_dir/ens_known_genes_DCPS_TIRAP_RP11-712L6.5.bed
#gene_bed=$gene_bed_total
run_reann_cluster $outdir $cff $gene_bed

# TEST FUSIONS
#HES1--IL29
cff=$reann_test_dir/HES1--IL29/HES1--IL29.cff
gene_bed=$reann_test_dir/ens_known_genes_HES1_IL29.bed
run_reann_cluster $outdir $cff $gene_bed
#SAA1--ATP1A1
cff=$reann_test_dir/SAA1--ATP1A1/SAA1--ATP1A1.cff
gene_bed=$reann_test_dir/ens_known_genes_SAA1--ATP1A1.bed
run_reann_cluster $outdir $cff $gene_bed
#BCL3--CTB-171A8.1
cff=$cff_dir/BCL3--CTB-171A8.1.cff
gene_bed=$bed_dir/BCL3--CTB-171A8.1.bed
run_reann_cluster $outdir $cff $gene_bed
#ANKIB1--AKAP9
cff=$reann_test_dir/cff_files/merged.cff_AKAP9--ANKIB1
gene_bed=$reann_test_dir/ens_known_genes.AKAP9--ANKIB1.bed
run_reann_cluster $outdir $cff $gene_bed
#FOSB--AADACL2
cff=$reann_test_dir/FOSB--AADACL2/FOSB--AADACL2.cff
gene_bed=$bed_dir/ens_known_genes.MIR548H2.FOSB.AADACL2.bed
run_reann_cluster $outdir $cff $gene_bed

exit 0

#SAA1--ATP1A1
outdir=$reann_test_dir/SAA1--ATP1A1
cff=$reann_test_dir/SAA1--ATP1A1/SAA1--ATP1A1.cff
#run_reann_cluster $outdir $cff $gene_bed_total

#SIM45.SIM52
cff=/hpf/largeprojects/ccmbio/mapostolides/validate_fusion/montreal_benchmark_validation/sim45.sim52.combined/output-7CALLERS-April-22-2020/fusions/cff/merged.cff
#outdir=/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/TEST_bwafilter-May-4-2020/SIM45.SIM52.cff_10000_bp.May-6-2020
outdir=$reann_test_dir/SIM45.SIM52.May-12-2020
run_reann_cluster $outdir $cff $gene_bed_total

#BEERS -ive control
cff=/hpf/largeprojects/ccmbio/mapostolides/validate_fusion/montreal_benchmark_validation/negative_control_beers/output-7CALLERS-April-17-2020/fusions/cff/merged.cff
outdir=/hpf/largeprojects/ccmbio/mapostolides/MODULES/RUN_BENCHMARKING_TOOLKIT/neg.control.beers.May-5-2020
#run_reann_cluster $outdir $cff $gene_bed_total

##BT474.KPL4.MCF7.SKBR3
outdir=/hpf/largeprojects/ccmbio/mapostolides/MODULES/RUN_BENCHMARKING_TOOLKIT/BT474.KPL4.MCF7.SKBR3.cff-April-30-2020
cff=/hpf/largeprojects/ccmbio/mapostolides/validate_fusion/haas_2019_validation/output.BT474.KPL4.MCF7.SKBR3-April-9-2020/fusions/cff/merged.cff
#run_reann_cluster $outdir $cff $gene_bed_total

##SIM50
outdir=/hpf/largeprojects/ccmbio/mapostolides/MODULES/RUN_BENCHMARKING_TOOLKIT/SIM50.cff-April-30-2020
cff=/hpf/largeprojects/ccmbio/mapostolides/validate_fusion/haas_2019_validation/output-SIM-April-17-2020/fusions/cff/merged.cff
#run_reann_cluster $outdir $cff $gene_bed_total

##
cff=$reann_test_dir/cff_files/merged.cff_BCL3 

#BCL3--CTB-171A8.1
cff=$cff_dir/BCL3--CTB-171A8.1.cff
gene_bed=$bed_dir/BCL3--CTB-171A8.1.bed
gene_bed=$gene_bed_total
#ANKIB1--AKAP9
#gene_bed=$reann_test_dir/ens_known_genes.ANKIB1--AKAP9.bed
#cff=$reann_test_dir/cff_files/merged.cff_AKAP9--ANKIB1
#run_reann_cluster $outdir $cff $gene_bed

#FOSB--AADACL2
cff=$reann_test_dir/FOSB--AADACL2/FOSB--AADACL2.cff
gene_bed=$bed_dir/ens_known_genes.MIR548H2.FOSB.AADACL2.bed
#1 entry, failing due to breakpt being 1 bp off:
#cff=$reann_test_dir/FOSB--AADACL2/FOSB--AADACL2.cff.test
#run_reann_cluster $outdir $cff $gene_bed

