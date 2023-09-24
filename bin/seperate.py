#!/usr/bin/env python3
# -- coding: utf-8 --
import sys
import os
infile=sys.argv[1]
outdir=sys.argv[2]
if not os.path.exists(outdir+'/ref'):
    os.mkdir(outdir+'/ref')

file=open(infile,'r')
for line in file:
    if line[0] == '>':
        name = line.strip().split(' ')[0]
        virus=name.split('>')[1]
        if '.' in virus: 
            print(virus)
            name1=virus.split(' ')[0]
            fw= open(outdir+'/ref/'+virus+'.fa','w')
            fw.write('>'+name1+'\n')
        else:
            
            fw= open(outdir+'/ref/'+virus+'.fa','w')
            fw.write('>'+virus+'\n')
    else:
            fw.write(line)
file.close()
