#!/usr/bin/env python3
# -- coding: utf-8 --
"""
python result.py $meta $rag_fra $rag_trinity_fra $trinity
@author: fp
Determine which genome graction is larger
"""
import pandas as pd
import numpy as np
import sys
import os
metacompass=sys.argv[1]
rag_fra=sys.argv[2]
rag_trinity_fra=sys.argv[3]
trinity=sys.argv[4]
totalreads=sys.argv[5]
virus=sys.argv[6]
sample=sys.argv[7]
if os.path.exists(metacompass):
    metacompass=open(metacompass,'r')
    next(metacompass)
    for line in metacompass:
        meta=line.strip().split(' ')[-1]
else:
    meta=str(0)

if os.path.exists(rag_fra):
    rag_f=open(rag_fra,'r')
    next(rag_f)
    for line1 in rag_f:
        rag=line1.strip().split(' ')[-1]
else:
    rag=str(0)

if os.path.exists(rag_trinity_fra):
    rag_trinity_file=open(rag_trinity_fra,'r')
    next(rag_trinity_file)
    for line2 in rag_trinity_file:
        rag_trinity=line2.strip().split(' ')[-1]
else:
    rag_trinity=str(0)

if os.path.exists(trinity):
    trinity=open(trinity,'r')
    next(trinity)
    for line3 in trinity:
        trinity=line3.strip().split(' ')[-1]
else:
    trinity=str(0)

dir=rag_fra.split('genome')[0]
out=rag_fra.split('metaquast')[0]
outfile=out+'result/'+virus
print(outfile)

if not os.path.exists(outfile):
     os.makedirs(outfile) #创建多级目录用makedirs，单击目录mkdir
resultfile=open(outfile+'/result.txt','w')
resultfile.write('ragtag_trinity(%)\t'+'ragtag(%)\t'+'metacompass(%)\t'+'Trinity(%)\t'+'abundance\t'+'coverage(%)\t'+'depth_cov\t'+'depth_all'+'\n')

if max(round(float(rag),3),round(float(rag_trinity),3),round(float(meta),3),round(float(trinity),3))==round(float(rag_trinity),3):
    print('Genome fraction: ragtag_trinity '+rag_trinity+' is the largest, and move to the result dir')
    resultfile.write(rag_trinity+'\t'+rag+'\t'+meta+'\t'+trinity+'\t')
    os.system('cat '+out+'ragtag_output_trinity/ragtag.scaffold.fasta1|grep '+virus+' -A 1 >'+outfile+'/result.fasta')
elif max(round(float(rag),3),round(float(rag_trinity),3),round(float(meta),3),round(float(trinity),3))==round(float(rag),3):
    print('Genome fraction: ragtag '+rag+' is  the largest , and move to the result dir')
    resultfile.write(rag_trinity+'\t'+rag+'\t'+meta+'\t'+trinity+'\t')
    os.system('cat  '+out+'ragtag_output/'+virus+'/ragtag.scaffold.fasta1|grep '+virus+' -A 1 > '+outfile+'/result.fasta')
elif max(round(float(rag),3),round(float(rag_trinity),3),round(float(meta),3),round(float(trinity),3))==round(float(meta),3):
    print('Genome fraction: metacompass '+meta+' is the largest, and move to the result dir')
    resultfile.write(rag_trinity+'\t'+rag+'\t'+meta+'\t'+trinity+'\t')
    os.system('cat '+out+'metacompass/metacompass_output/metacompass.final.ctg.fa1|grep '+virus+' -A 1 > '+outfile+'/result.fasta')
elif max(round(float(rag),3),round(float(rag_trinity),3),round(float(meta),3),round(float(trinity),3))==round(float(trinity),3):
    print('Genome fraction: Trinity '+trinity+' is the largest, and move to the result dir')
    resultfile.write(rag_trinity+'\t'+rag+'\t'+meta+'\t'+trinity+'\t')
    os.system('cp '+dir+'blastvirus_fa/' +sample+'.fa '+outfile+'/result.fasta')

##abundance

virusbam = os.path.getsize(out+sample+'.sorted.bam') 
if virusbam<100:
    print('virusbam size:'+str(virusbam)+'. It may be an empty file')
else:
    print('virusbam size:'+str(virusbam))

maxnum=max(round(float(rag),3),round(float(rag_trinity),3),round(float(meta),3),round(float(trinity),3))


if int(virusbam) !=int(62):
    bamfile=os.popen('samtools idxstats '+out+sample+".sorted.bam|awk 'NR>1{print line}{line=$0}'")#删除尾行
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
                print('totalreads:\t'+totalreads+'\n'+'viruslen:\t',len,'\n','virus true len:\t',lentrue,'\n','mappedreads:\t',mappedreads,'\n','abundance:\t',abundance,'\n')
                resultfile.write(str(abundance))
else:
    print(virusbam,' maybe empty')


os.system("bedtools genomecov -ibam "+out+sample+'.sorted.bam'+" -bga > "+out+"result/"+"genomecov.txt")

inputfile=out+"result/"+"genomecov.txt"

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

