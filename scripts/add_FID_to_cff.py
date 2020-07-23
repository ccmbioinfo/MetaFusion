#!/usr/bin/env python
import sys
#sys.path.append("/hpf/largeprojects/ccmbio/jiangyue/DIPG_analysis_by_samples/Scripts/pygeneann/pygenefusionann")
sys.path.append("/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin")
import  pygeneann_reads_capture_DEV as pygeneann
import sequtils
import pysam
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('cff_file', action='store', help='CFF file, can be .cff or cff.reann')
args = parser.parse_args()
cff_file = args.cff_file

#cff_file="/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/testing_pipeline/UHRR.benchmark.June-15-2020/tmp"
n = 1
for line in open(cff_file, "r"):
    fusion = pygeneann.CffFusion(line)

    #annotate fusion id and seq
    fusion.fusion_id = "F" + (str(n)).zfill(8)
    print fusion.tostring()
    n += 1


