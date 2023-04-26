#!/usr/bin/env python
"""
@author: fuping
source /home/root640/software/miniconda3/bin/activate fp
"""

import os
import argparse
import time
import sys
parser = argparse.ArgumentParser()
parser.add_argument('--evalue', help='The BLAST E-value. (default: 1e-5)', default='1e-5')
requiredNamed = parser.add_argument_group('required arguments')
requiredNamed.add_argument('--fastq_1', help='raw fastq.gz file paired-end 1', required=True)
requiredNamed.add_argument('--fastq_2', help='raw fastq.gz file paired-end  2', required=True)
#requiredNamed.add_argument('--fastq', help='raw fastqgz file with single-end', required=True)
requiredNamed.add_argument('--outdir', help='output path', required=True)
requiredNamed.add_argument('--Diamond_VirusProtein_db', help='RefSeqVirusProtein db maked by diamond', required=True)
requiredNamed.add_argument('--Diamond_nr_db', help='raw fastqgz file with single-end', required=True)
requiredNamed.add_argument('--threads', help='running threads (default: 20)', default=20, type=int)
requiredNamed.add_argument('--clean_fastq_1', help='clean fastqgz file paired 1', required=True)
requiredNamed.add_argument('--clean_fastq_2', help='clean fastqgz file paired 2', required=True)
requiredNamed.add_argument('--MetaCompass_dir', help='MetaCompass software path', required=True)
args = parser.parse_args()

softdir = sys.path[0]
evalue=str(args.evalue)
fastq_1 = args.fastq_1.split('.gz')[0]
fastq_2 = args.fastq_2.split('.gz')[0]
outdir = args.outdir
Diamond_VirusProtein_db= args.Diamond_VirusProtein_db
Diamond_nr_db= args.Diamond_nr_db
threads= str(args.threads)
virus_fasta= args.virus_fasta
MetaCompass_dir=args.MetaCompass_dir
name1= args.fastq_1.split('/')[-1]
name= name1.split('_1.fastq.gz')[0]
clean_fastq_1=args.clean_fastq_1
clean_fastq_2=args.clean_fastq_2
Trinity_out=outdir+"/Trinity/trinity_out_"+name
Diamond_out=outdir+"/Diamond/"+name
Classify_out=outdir+"/Classify/"+name
if not os.path.exists(outdir):
    os.makedirs(outdir) 
if not os.path.exists(outdir+"/Fastp/"):
    os.makedirs(outdir+"/Fastp/") 
if not os.path.exists(outdir+"/Diamond/"):
    os.makedirs(outdir+"/Diamond/") 
if not os.path.exists(Trinity_out):
    os.makedirs(Trinity_out) 
if not os.path.exists(outdir+"/Classify/"):
    os.makedirs(outdir+"/Classify/") 
if not os.path.exists(fastq_1 and fastq_2):
    os.system('gunzip '+args.fastq_1+";"+'gunzip '+args.fastq_2)

if not os.path.exists(softdir+''/../db/final_out.fa'):
    os.system('gunzip '+softdir+''/../db/final_out.fa.gz')
    
if not os.path.exists(softdir+'/../db/virus.genomic.fna.nhr'):
    os.system('makeblastdb  -dbtype nucl  -in '+softdir+'/../db/final_out.fa  -input_type fasta  -parse_seqids  -out  '+softdir+'/../db/virus.genomic.fna')

if not os.path.exists(outdir+'/Ref/'+name+'.fa'):
    os.system('python '+softdir+'/0_run1_paired.py  --fastq_1 '+fastq_1+' --fastq_2 '+fastq_2+' --outdir '+outdir+' --Diamond_VirusProtein_db '+Diamond_VirusProtein_db+' --Diamond_nr_db '+Diamond_nr_db+' --evalue '+evalue+' --threads '+threads)
    os.system('gzip -k -d '+clean_fq1+';gzip -k -d '+clean_fq2)
    os.system('python '+softdir+'/0_run2.py '+' --clean_fastq_1 '+clean_fastq_1+'.gz   --clean_fastq_2 '+clean_fastq_2+'.gz --virus_fasta '+outdir +'/Ref/'+name+'.fa --outdir '+outdir+' --MetaCompass_dir '+MetaCompass_dir+' --threads 20')
if  os.path.exists(outdir+'/Ref/'+name+'.fa'):
    os.system('python '+softdir+'/0_run2.py '+' --clean_fastq_1 '+clean_fastq_1+'  '+' --clean_fastq_2 '+clean_fastq_2+' --virus_fasta '+virus_fasta+' --outdir '+outdir+' --MetaCompass_dir '+MetaCompass_dir+' --threads 20')
