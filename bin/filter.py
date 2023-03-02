#!/usr/bin/env python3
# -- coding: utf-8 --

"""
@author: fp
#2021.6.11
"""
import pandas as pd
import numpy as np
import sys
import os

besthit_doc=sys.argv[1]
virussure_doc=sys.argv[2]
species_doc=sys.argv[3]
genus_doc=sys.argv[4]
family_doc=sys.argv[5]
name=sys.argv[6]
software=sys.argv[7]
def judgmentvirus():
#read all diamond blast besthit
    
    dfall=pd.read_csv(besthit_doc,sep='\t',header=None, names=['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore','slen','stitle','salltitles','qcovhsp','nident','staxids'])
    dfall1=dfall.dropna(subset=['staxids'])
#read all virus taxid
    dfvirus=pd.read_csv(software+'/db/virus20220308.all.prot2tax',sep='\t',header=None, names=['acc','accversion','taxid','GI'])
    viruslist=dfvirus['taxid'].drop_duplicates().values.tolist()
    viruslist1= list(map(str,viruslist))
#get overlap
    dfsure=dfall1[dfall1.staxids.isin(viruslist1)]
    dfsure.to_csv(virussure_doc, index=False, header=None,sep='\t')

def filter():        
        df1=pd.read_csv(virussure_doc,sep='\t',header=None, names=['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore','slen','stitle','salltitles','qcovhsp','nident','staxids'])
        #print(df1)#      staxids qcovhsp
        df1['staxids'] = df1['staxids'].astype(int)
        df2 = df1[(df1['pident']>=90)&(df1['qcovhsp']>=80)]
        df2.to_csv(species_doc, index=False, header=None,sep='\t')
        df3 = df1[(df1['pident']>=70)&(df1['qcovhsp']>=60)]
        genus = pd.concat([df3, df2, df2]).drop_duplicates(keep=False)
        genus.to_csv(genus_doc, index=False, header=None,sep='\t')
        df4 = df1[(df1['pident']>=50)&(df1['qcovhsp']>=40)]
        family=pd.concat([df4, df2, df2]).drop_duplicates(keep=False)
        family1=pd.concat([df4, df3, df3]).drop_duplicates(keep=False)
        family1.to_csv(family_doc, index=False, header=None,sep='\t')

def getrefseq():
        dic1={}
        taxidset=[]
        notin=[]
        file=open(species_doc,'r')
        infile=open(software+'/db/taxid_ref_3','r')
        final=open(species_doc+'final.txt','w')
        for line in infile:
            taxid=line.strip().split('\t')[0]
            refseqid=line.strip().split('\t')[1]
            dic1[taxid]=[refseqid]
        for i in file:
            taxid2=i.strip().split('\t')[-1]
            taxid1=taxid2.rstrip('')
            taxidset.append(taxid1)
        taxidset=sorted(list(set(taxidset)))
        print('########################################################################\n')
        print('The identified virus id is as follows:\n',taxidset,'\nnumbers of virus species: ',len(taxidset),'\n')
        print('########################################################################')
        for each in taxidset:
            if each in dic1.keys():
                 print(each+' refseqid is '+str(list(set(dic1[each]))))
                 lin='\t'.join(os.popen("echo "+each+" |taxonkit lineage |taxonkit reformat  -P |awk -F '\t' '{print $3}'"))
                 print(each+'\tlineage:\t'+str(lin))
                 final.write('\t'.join(list(set(dic1[each])))+'\t'+each+'\t'+str(lin))
            else:
                notin.append(each)
        #print('The virus not in our database :',len(notin))
        if len(notin)==0:
            pass
        else:
            print(', '.join(notin)+' is not in our database, please download the sequence from NCBI and put it in the ./genome/virus/ directory and named as '+name+'.fa')
 

def getfa():
      refseq=open(species_doc+'final.txt','r')
      if os.path.exists('./genome/virus/' +name +'.fa'):
          os.remove('./genome/virus/' +name +'.fa')
      for line in refseq:
          ref=line.strip().split('\t')[0] 
          #taxid=line.strip().split('\t') [1] 

          if ', ' in ref:
              ref1=ref.split(', ')
              for i in ref1:
                  print('########################################################################')
                  print(i,'    get fasta to ./genome/virus/ directory')
                  #print('contain two or more refgenome:')
                  os.system("blastdbcmd -db "+software+"/db/virushostdb.genomic.fna  -dbtype nucl -entry "+i+' > ./genome/virus/' +i +'.fa') 
                  os.system("cat ./genome/virus/"+i+'.fa >> ./genome/virus/' +name +'.fa') 
          else:
              os.system("blastdbcmd -db "+software+"/db/virushostdb.genomic.fna  -dbtype nucl -entry "+ref+' > ./genome/virus/' + ref+'.fa')
              os.system("cat ./genome/virus/"+ref+'.fa >> ./genome/virus/' +name +'.fa') 
              print('########################################################################')
              print(ref,'    get fasta to ./genome/virus/ directory '+name+'.fa')

if  __name__=='__main__':
    judgmentvirus()
    filter()
    getrefseq()
    getfa()
