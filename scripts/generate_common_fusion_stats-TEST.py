#!/usr/bin/env python
import sys
#sys.path.append("/hpf/largeprojects/ccmbio/jiangyue/DIPG_analysis_by_samples/Scripts/pygeneann/pygenefusionann")
#sys.path.append("/hpf/largeprojects/ccmbio/jiangyue/Genap_ccm/pygenefusionann/")
import pygeneann
import sequtils
import argparse
import os
#parser = argparse.ArgumentParser()
#parser.add_argument('cff_file', action='store', help='CFF file')

#args = parser.parse_args()

top_dir="/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/"
#cff_file="/Users/mapostolides/Desktop/mugqic_tools/python-tools/fusiontools/0.1.0/bin/reann_cff_fusion_testing/cff_files/merged.cff_RPRD2--LAMC2.reann"
#cff_file="/Users/mapostolides/Desktop/mugqic_tools/python-tools/fusiontools/0.1.0/bin/reann_cff_fusion_testing/cff_files/merged.cff_RPRD2--LAMC2-6fusions.reann"
cff_file=os.path.join(top_dir, "reann_cff_fusion_testing/cff_files/merged.cff.reann")

cffstats = pygeneann.CffFusionStats(cff_file)	

cffstats.generate_common_fusion_stats_by_genes(cff_file)

