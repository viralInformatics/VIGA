#!/usr/bin/env python3
# -- coding: utf-8 --
"""

@author: fp
Determine which genome graction is larger
"""
import pandas as pd
import numpy as np
import sys
import os

virus=sys.argv[1]
rag_fra=sys.argv[2]
totalreads=sys.argv[3]
outdir=sys.argv[4]


if os.path.exists(rag_fra):
    rag_f=open(rag_fra,'r')
    next(rag_f)
    for line1 in rag_f:
        rag=line1.strip().split(' ')[-1]
else:
    rag=str(0)

virus=str(virus)
outfile=outdir+'/result/'+virus
print(outfile)

if not os.path.exists(outfile):
     os.makedirs(outfile) #创建多级目录用makedirs，单击目录mkdir
resultfile=open(outfile+'/result.txt','w')
resultfile.write('VIGA(%)\t'+'abundance\t'+'coverage(%)\t'+'depth_cov\t'+'depth_all'+'\n')

if round(float(rag),3)!=0:
    print('************\nVIGA Genome fraction: '+rag+', and move to the result dir')
    resultfile.write(rag+'\t')
    os.system('cat  '+outdir+'/ragtag_output/'+virus+"/ragtag.scaffold.fasta1|grep "+virus+" -A 1 > "+outfile+"/result.fasta")

##abundance

virusbam = os.path.getsize(outdir+'/virus_sorted.bam') 
if virusbam<100:
    print('virusbam size:'+str(virusbam)+'. It may be an empty file')
else:
    print('virusbam size:'+str(virusbam))

maxnum=round(float(rag),3)

if int(virusbam) !=int(62):
    os.system("samtools index "+outdir+'/virus_sorted.bam')
    bamfile=os.popen('samtools idxstats '+outdir+"/virus_sorted.bam|awk 'NR>1{print line}{line=$0}'")#删除尾行
    for line in bamfile:
        if virus in line:
            len=line.strip().split('\t')[1]
            lentrue=int(len)*maxnum*0.01
            mappedreads=line.strip().split('\t')[2]
            di=float(int(totalreads)*lentrue)
            if di==0:
                abundance=0.0001
                resultfile.write(str(abundance))
            else:
                abundance=round(float(mappedreads)*1000*1000000/di,3)
                #print('totalreads:\t'+totalreads+'\n'+'viruslen:\t',len,'\n','virus true len:\t',lentrue,'\n','mappedreads:\t',mappedreads,'\n','abundance:\t',abundance,'\n')
                resultfile.write(str(abundance))
else:
    print(virusbam,' maybe empty')


os.system("bedtools genomecov -ibam "+outdir+'/virus_sorted.bam'+" -bga > "+outdir+"/result/"+"genomecov.txt")

inputfile=outdir+"/result/"+"genomecov.txt"

inputFile = open(inputfile, "r")
length = 0
covLength = 0
counts = 0
for line in inputFile:
	lineAS = line.split("\t")
	chromosome = lineAS[0]
	start = lineAS[1]
	end = lineAS[2]
	count = lineAS[3].split("\n")[0]
	if virus in chromosome:
	    length_tmp = int(end) - int(start)
	    length += length_tmp
	    counts += (int(float(count)) * length_tmp)
	    if count != "0":
                            covLength_tmp = length_tmp
                            covLength += covLength_tmp
if float(length)>0 and float(covLength)>0:
    coverage = "%.2f" % ((float(covLength) / float(length)) * 100)
    depth_cov = "%.2f" % (float(counts) / float(covLength))
    depth_all = "%.2f" % (float(counts) / float(length))
else:
    coverage = 0
    depth_cov = 0
    depth_all = 0


inputFile.close()
resultfile.write('\t'+str(coverage)+'\t'+str(depth_cov)+'\t'+str(depth_all)+'\n')
resultfile.close()
os.system('rm '+outdir+"/result/"+"genomecov.txt")

