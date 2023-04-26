#!/usr/bin/env python3
# -- coding: utf-8 --

"""
@author: fp
"""

import pandas as pd
import numpy as np
import sys
import os

srr=sys.argv[1]
software=sys.argv[2]
outdir=sys.argv[3]
genus=sys.argv[4]

dic={}
alltaxid=[]
allgenusid=[]

path=outdir+'/newvirus/genus_result/'
index=outdir+'/newvirus/index/'
if not os.path.exists(path):
    os.makedirs(path) 
if not os.path.exists(index):
    os.makedirs(index) 

def getgenus():
    global genus
    global family

    df=pd.read_csv(genus,sep='\t',header=None, names=['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore','slen','stitle','salltitles','qcovhsp','nident','staxids'])
    df['staxids']=df['staxids'].astype("int")#.map(lambda x: str(x)[:-2])#delete the last 2
    df['staxids']=df['staxids'].astype("str")
    tmp={k: v['qseqid'].tolist() for k, v in df.groupby("staxids")}
    for each in tmp.keys():
        each=str(each)
        lin='\t'.join(os.popen("echo "+each+" |taxonkit lineage |taxonkit reformat  -P |awk -F '\t' '{print $3}'"))
        genus=lin.strip().split(';')[-2]
        family=lin.strip().split(';')[-3]
        tmp[each].append('|'+genus)
        tmp[each].append('|'+family)
    print(tmp)

    getfa=pd.DataFrame(list(tmp.items()),columns=['taxid','contig|genus|family'])
    pd.set_option('display.max_rows', None)
    getfa['contig']=getfa['contig|genus|family'].astype('str').map(lambda x:x.split("', '|")[0])
    getfa['contig']=getfa.contig.str.replace('[','',regex=True)
    getfa['contig']=getfa.contig.str.replace("'",'',regex=True)

    getfa['genus']=getfa['contig|genus|family'].astype('str').map(lambda x:x.split("', '|")[-2])
    #getfa['genus']=getfa['genus'].replace("', '",'',inplace = True)
    getfa['family']=getfa['contig|genus|family'].astype('str').map(lambda x:x.split("', '|")[-1])
    getfa.family=getfa.family.str.replace("']",'',regex=True)
    getfa1=getfa[['genus','family','contig']]
    getfa1.to_csv(outdir+'/'+srr+'_tmp',index=False,header=0,sep='\t')

def getfa():
    infile=open(outdir+'/'+srr+'_tmp','r')
    global genusid
    if not os.path.exists(path):
        os.makedirs(path)
    if not os.path.exists(index):
        os.makedirs(index)
    for line in infile:
        genusid=line.strip().split('\t')[0]
        contig=line.split('\t')[-1]
        contig1=contig.split(', ')
        if os.path.exists(path+srr+'_'+genusid+'_anno.txt'):
            os.remove(path+srr+'_'+genusid+'_anno.txt')
            os.remove(path+srr+'_'+genusid+'.fa')
        os.system("echo 'contig|protein|pident|length|protein title|qcovhsp' >> " +path+srr+'_'+genusid+'_anno.txt' )
        os.system('seqkit seq -w 0 '+outdir+'/Trinity/trinity_out_'+srr+'/Trinity.fasta > '+outdir+'/Trinity/trinity_out_' +srr+'/Trinity.fasta1'+ "> /dev/null 2>&1")
        for i in contig1:
            i=i.strip()
            os.system('grep -A 1 '+i+'  '+outdir+'/Trinity/trinity_out_' +srr+'/Trinity.fasta1 >>'+path+srr+'_'+genusid+'.fa1')
            os.system('cat '+path+srr+'_'+genusid+'.fa1 |'+' seqkit rmdup -s -o '+path+srr+'_'+genusid+'.fa')
            
            os.system('grep '+i+' '+outdir+'/Classify/'+srr+'.genus'+" |awk -F '\t'  -v OFS='|'  '{print $1,$2,$3,$4,$14,$16}'"+' >> '+path+srr+'_'+genusid+'_anno.txt')
        allgenusid.append(genusid)
        
    os.remove(outdir+'/Trinity/trinity_out_' +srr+'/Trinity.fasta1')
    allgenusid_=list(set(allgenusid))
    return allgenusid_

def getcoverage():
    allgenusid_u=getfa()
    print(allgenusid_u)
    genusdic={}
    genus_len=open(software+'/../db/genus_len','r')
    for line in genus_len:
        genus1=line.strip().split('\t')[0]
        len=line.strip().split('\t')[1]
        genusdic[genus1]=[len]
    for oneofgenus in allgenusid_u:
        if oneofgenus in genusdic.keys():
            genusfa=path+srr+'_'+oneofgenus+'.fa'
            num=float("".join([str(x) for x in genusdic[oneofgenus]])) #float() argument must be a string or a number, not 'list'
            avglen=round(num,2)
        
            os.system('bwa index -a bwtsw -p '+index+srr+'_'+oneofgenus+' '+genusfa)
            os.system('bwa mem -T 20 -t 20 '+index+srr+'_'+oneofgenus+ ' '+outdir+'/Fastp/'+srr+'_1.clean.fastq ' + outdir+'/Fastp/'+srr+'_2.clean.fastq   | samtools sort |samtools view -F4 -o '+ index+srr+'_'+oneofgenus+'.sorted.bam')
            os.system('samtools index -@ 20 '+index+srr+'_'+oneofgenus+'.sorted.bam')
            os.system("bedtools genomecov -ibam "+index+srr+'_'+oneofgenus+'.sorted.bam'+" -bga > "+path+srr+'_'+oneofgenus+".genomecov.txt")

            inputfile=path+srr+'_'+oneofgenus+".genomecov.txt"
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
                length_tmp = int(end) - int(start)
                length += length_tmp
                counts += (int(float(count)) * length_tmp)
                if count != "0":
                    covLength_tmp = length_tmp
                    covLength += covLength_tmp
            if float(length)>0 and float(covLength)>0:
                coverage = "%.2f" % ((float(covLength) / float(length)) * 100)
            else:
                coverage = 0
            result=path+srr+'_'+oneofgenus+"_anno_detail.txt"
            resultfile=open(result,'w')
            resultfile.write('Reads Coverage(%)\tGenome Fraction(%)\tAbudance(RVKM)\n')
       ########Reads Coverage(%)
            resultfile.write(str(coverage)+'\t')
       ########Genome Fraction(%)
            contigsum=os.popen('samtools idxstats '+index+srr+'_'+oneofgenus+".sorted.bam|awk '{print $2}' | awk '{sum+=$1}END{print sum}'")#samtools idxstats第二列是序列长度的加和
            contigsumlen=contigsum.read()
 
            genusfraction="%.2f" % ((int(contigsumlen) / float(avglen)) * 100)
            resultfile.write(str(genusfraction)+'\t')
        ########Abudance(RVKM)
            mapped=os.popen('samtools idxstats '+index+srr+'_'+oneofgenus+".sorted.bam|awk '{print $3}' | awk '{sum+=$1}END{print sum}'")#samtools idxstats第三列是mapped reads数
            mappedreads=mapped.read()
            total=os.popen('cat '+outdir+'/Fastp/'+srr+'_1.clean.fastq |wc -l')
            totalreads=total.read()
            print('totalreads\t','contigsumlen\t','avglen\t')
            print(int(totalreads)/2,int(contigsumlen),float(avglen))
            di=float((int(totalreads)/2)*int(contigsumlen)*float(genusfraction))
            if di==0:
                abundance=0
            else:
                abundance=round(float(mappedreads)*1000*10000/di,3)
#分子去掉了100，因为genusfraction那里是百分数形式
        ###########result
            resultfile.write(str(abundance)+'\n')
            os.remove(index+srr+'_'+oneofgenus+'.sorted.bam')
            os.remove(index+srr+'_'+oneofgenus+'.sorted.bam.bai')
           #os.remove(path+srr+'_'+oneofgenus+".genomecov.txt")
    os.remove(outdir+'/'+srr+'_tmp')
    os.system('rm '+path+'/*fa1')

if  __name__=='__main__':
    getgenus()
    getfa()
    getcoverage()
