#!/usr/bin/env python
import sys
import pygeneann
import sequtils
import pysam
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('cff_file', action='store', help='CFF file, can be .cff or cff.reann')
parser.add_argument('bwa_sam', action='store', help='bwa alignment of the reann file converted fasta')

args = parser.parse_args()

samf = open(args.bwa_sam, "r")
f = pysam.AlignmentFile(samf, 'r')

for line in open(args.cff_file, "r"):
    fusion = pygeneann.CffFusion(line)
    n = 0
    score = 0
    read1 = ""
    read2 = ""
    for read in f.fetch(until_eof=True):
        if not read1:
            read1 = read
        elif not read2:
            read2 = read
            fusion_id1 = read1.query_name.split("_")[0]    
            fusion_id2 = read2.query_name.split("_")[0]
            assert fusion_id1 == fusion_id2 == fusion.fusion_id, fusion_id1 + fusion_id2 + fusion.fusion_id
            if not (read1.is_unmapped or read2.is_unmapped):
                as_1 = int(read1.get_tag("AS"))
                as_2 = int(read2.get_tag("AS"))
                xs_1 = int(read1.get_tag("XS"))
                xs_2 = int(read2.get_tag("XS"))
                if as_1 > xs_1 and as_2 > xs_2:
                    print fusion.tostring()
                
            read1 = ""
            read2 = ""
            break
                
            
    '''
    while True:
        line2 = samf.readline()
        if line2.startswith("@"):
            continue
        n += 1
        tmp = line2.split()
        fusion_id = tmp[0].split("_")[0]
        if fusion_id != fusion.fusion_id:
            print >> sys.stderr, "cff and sam not match."
            sys.exit(1)
        score += int(tmp[4])
        if (n == 2):
            break    
    if score == 120:
        print fusion.tostring()
    '''
samf.close()    
f.close()
