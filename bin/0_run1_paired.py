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
args = parser.parse_args()

softdir = sys.path[0]
evalue=str(args.evalue)
fastq_1 = args.fastq_1.split('.gz')[0]
fastq_2 = args.fastq_2.split('.gz')[0]
outdir = args.outdir
Diamond_VirusProtein_db= args.Diamond_VirusProtein_db
Diamond_nr_db= args.Diamond_nr_db
threads= str(args.threads)

name1= args.fastq_1.split('/')[-1]
name= name1.split('_1.fastq.gz')[0]
clean_out=outdir+"/Fastp/"+name
clean_fq1=clean_out+'_1.clean.fastq'
clean_fq2=clean_out+'_2.clean.fastq'
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

if os.path.exists(softdir+'/../db/virus.genomic.fna.nhr'):
    os.system('makeblastdb  -dbtype nucl  -in '+softdir+'/../db/final_out.fa  -input_type fasta  -parse_seqids  -out  '+softdir+'/../db/virus.genomic.fna')



if os.path.exists(fastq_1 and fastq_2):
    print('*********fastp begins :', time.asctime(time.localtime(time.time())),'\n')
    os.system('fastp -i '+fastq_1+' -I '+fastq_2+'  -o '+clean_fq1+' -O '+clean_fq2+"  -h "+clean_out+'.html -w '+threads+' -j '+clean_out+'.json')
    os.system("sed -i 's/ //g' "+clean_fq1)
    os.system("sed -i 's/ //g' "+clean_fq2)
else:
    print('*********fastp failed :', time.asctime(time.localtime(time.time())),'\n')
if os.path.exists(clean_fq1):
    print('*********Trinity denovo assembly begins :', time.asctime(time.localtime(time.time())),'\n')
    os.system("Trinity --seqType fq --left "+clean_fq1+"  --right "+clean_fq2+" --output "+Trinity_out+" --CPU 26 --max_memory 50G")
else:
    print('*********Trinity failed :', time.asctime(time.localtime(time.time())),'\n')
if os.path.exists(Trinity_out+'/Trinity.fasta'):
    print('*********diamond blastx against virus protein database begins :', time.asctime(time.localtime(time.time())),'\n')
    os.system("diamond blastx -q "+Trinity_out+'/Trinity.fasta'+ ' --db '+Diamond_VirusProtein_db+' -o '+Diamond_out+'.vp.txt -e '+evalue+'  --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore  slen  stitle salltitles qcovhsp nident staxids')
    os.system("cat "+Diamond_out+".vp.txt |awk '{print $1}'|uniq> "+Diamond_out+"vp1nd")
    os.system("cat "+Trinity_out+'/Trinity.fasta | seqkit grep -f  '+Diamond_out+"vp1nd>"+Diamond_out+".fa")
else:
    print('*********Trinity failed :', time.asctime(time.localtime(time.time())),'\n')
if os.path.exists(Diamond_out+".fa"):
    print('*********diamond blastx agaisnt nr begins :', time.asctime(time.localtime(time.time())),'\n')
    os.system("diamond blastx -q "+Diamond_out+".fa"+ ' --db '+Diamond_nr_db+' -o '+Diamond_out+'.nr.txt -e '+evalue+'  --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore  slen  stitle salltitles qcovhsp nident staxids')
    os.system("cat "+Diamond_out+'.nr.txt | sort -k1,1 -k12,12gr -k11,11g -k3,3gr | sort -u -k1,1 --merge>'+Diamond_out+'.besthit.txt')
    os.system("python  "+softdir+'/filter.py '+Diamond_out+'.besthit.txt '+Diamond_out+'.virussure.txt '+Classify_out+' '+name+' '+softdir+'/.. '+outdir)
else:
    print('*********Diamond failed :', time.asctime(time.localtime(time.time())),'\n')

if os.path.exists(Classify_out+".genus"):
    print('*********Novol virus discovery :', time.asctime(time.localtime(time.time())),'\n')
    os.system('python '+softdir+'/getLCA.py '+name+' '+softdir+' '+outdir+' '+Classify_out+".genus")
else:
    print(Classify_out+".genus")
    print('*********No genus fasta, Novol virus discovery finished :', time.asctime(time.localtime(time.time())),'\n')

if os.path.exists(Diamond_out+"vp1nd"):
    os.system('rm  '+Diamond_out+"vp1nd")

if not os.path.exists(clean_fq1+'.gz'):
    os.system('gzip '+clean_fq1)
    os.system('gzip '+clean_fq2)
print('*********Virus Identification and taxonomy finished', time.asctime(time.localtime(time.time())),'\n')



