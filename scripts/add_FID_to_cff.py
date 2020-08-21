#!/usr/bin/env python
import sys
import  pygeneann_MetaFusion as pygeneann
import sequtils
import pysam
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('cff_file', action='store', help='CFF file, can be .cff or cff.reann')
args = parser.parse_args()
cff_file = args.cff_file

n = 1
for line in open(cff_file, "r"):
    fusion = pygeneann.CffFusion(line)

    #annotate fusion id and seq
    fusion.fusion_id = "F" + (str(n)).zfill(8)
    print fusion.tostring()
    n += 1


