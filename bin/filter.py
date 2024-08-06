#!/usr/bin/env python3
# -- coding: utf-8 --

"""
@author: fp
"""
import pandas as pd
import numpy as np
import sys
import os

besthit_doc=sys.argv[1]
virussure_doc=sys.argv[2]
claassify_doc=sys.argv[3]
name=sys.argv[4]
software=sys.argv[5]
outdir=sys.argv[6]
def judgmentvirus():
#read all diamond blast besthit
    
    dfall=pd.read_csv(besthit_doc,sep='\t',header=None, names=['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore','slen','stitle','salltitles','qcovhsp','nident','staxids'])
    dfall1=dfall.dropna(subset=['staxids'])
    dfall1['staxids'] = dfall1['staxids'].astype(str)
#read all virus taxid
    dfvirus=pd.read_csv(software+'/db/virustaxid',sep='\t',header=None, names=['taxid'])
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
        df2.to_csv(claassify_doc+'.species', index=False, header=None,sep='\t')
        df3 = df1[(df1['pident']>=70)&(df1['qcovhsp']>=60)]
        genus = pd.concat([df3, df2, df2]).drop_duplicates(keep=False)
        genus.to_csv(claassify_doc+'.genus', index=False, header=None,sep='\t')
        df4 = df1[(df1['pident']>=50)&(df1['qcovhsp']>=40)]
        family=pd.concat([df4, df2, df2]).drop_duplicates(keep=False)
        family1=pd.concat([df4, df3, df3]).drop_duplicates(keep=False)
        family1.to_csv(claassify_doc+'.family', index=False, header=None,sep='\t')

def getrefseq():
        dic1={}
        taxidset=[]
        notin=[]
        file=open(claassify_doc+'.species','r')
        infile=open(software+'/db/taxid_ref_3','r')
        final=open(claassify_doc+'.final.txt','w')
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
      refseq=open(claassify_doc+'.final.txt','r')
      if not os.path.exists(outdir +'/Ref/'):
            os.makedirs(outdir +'/Ref/') 
      if os.path.exists(outdir +'/Ref/'+name +'.fa'):
          os.remove(outdir +'/Ref/'+name +'.fa')
      for line in refseq:
          ref=line.strip().split('\t')[0] 
          #taxid=line.strip().split('\t') [1] 

          if ', ' in ref:
              ref1=ref.split(', ')
              for i in ref1:
                  print('########################################################################')
                  print(i,' get fasta to '+outdir +'/Ref/'+name +' directory')
                  #print('contain two or more refgenome:')
                  os.system("blastdbcmd -db "+software+"/db/virus.genomic.fna  -dbtype nucl -entry "+i+' > '+outdir +'/Ref/' +i +'.fa') 
                  os.system('cat '+outdir +'/Ref/' +i +'.fa >> '+outdir+'/Ref/' +name +'.fa') 
          else:
              os.system("blastdbcmd -db "+software+"/db/virus.genomic.fna  -dbtype nucl -entry "+ref+' > '+outdir +'/Ref/'+ ref+'.fa')
              os.system('cat '+outdir +'/Ref/' +ref+'.fa >> '+outdir+'/Ref/' +name +'.fa') 
              print('########################################################################')
              print(ref,' get fasta to '+outdir +'/Ref/'+name +' directory')

if  __name__=='__main__':
    judgmentvirus()
    filter()
    getrefseq()
    getfa()
