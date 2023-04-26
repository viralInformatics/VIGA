#!/usr/bin/env python
"""
@author: fuping
source /home/root640/software/miniconda3/bin/activate fp
"""
import os
import argparse
import time
import sys
import pandas as pd
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('--len', help='length of reads (default: 100)', default='100')

requiredNamed = parser.add_argument_group('required arguments')
requiredNamed.add_argument('--clean_fastq_1', help='clean fastqgz file paired 1', required=True)
requiredNamed.add_argument('--clean_fastq_2', help='clean fastqgz file paired 2', required=True)
requiredNamed.add_argument('--virus_fasta', help='viral accession or chromosome, eg:NC_009334.1 or chr1,chr2, sure NC_009334.1,fa is in ./genome/virus/', required=True)
requiredNamed.add_argument('--outdir', help='output path', required=True)
requiredNamed.add_argument('--MetaCompass_dir', help='MetaCompass software path', required=True)
requiredNamed.add_argument('--threads', help='running threads (default: 20)', default=20, type=int)
args = parser.parse_args()


readslen=str(args.len)
clean_fastq_1 = args.clean_fastq_1.split('.gz')[0]
clean_fastq_2 = args.clean_fastq_2.split('.gz')[0]
virus = args.virus_fasta.split('.')[0]
virus_fasta= args.virus_fasta
outdir = args.outdir
MetaCompass_dir= args.MetaCompass_dir
virusfasta=args.virus_fasta
threads= str(args.threads)

softdir = sys.path[0]


virus_sorted_bam =outdir+'/virus_sorted.bam'
pipeup = outdir+'/virus.pipeup'
consensus = outdir+'/virus.consensus.fa'
metacompassout=outdir+'/metacompass'
ragtag_output=outdir+'/ragtag_output'
genome=outdir+'/metacompass'+'/metacompass_output/metacompass.final.ctg.fa'
genome1=outdir+'/metacompass'+'/metacompass_output/metacompass1.final.ctg.fa'
contig=outdir+'/metacompass'+'/assembly/contigs.fasta'
contig1=outdir+'/metacompass'+'/assembly/contigs1.fasta'

if not os.path.exists(softdir+''/../db/final_out.fa'):
    os.system('gunzip '+softdir+''/../db/final_out.fa.gz')
if not os.path.exists(softdir+'/../db/virus.genomic.fna.nhr'):
    os.system('makeblastdb  -dbtype nucl  -in '+softdir+'/../db/final_out.fa  -input_type fasta  -parse_seqids  -out  '+softdir+'/../db/virus.genomic.fna')

if not os.path.exists(outdir):
    os.mkdir(outdir)
if not os.path.exists(clean_fastq_1 and clean_fastq_2):
    os.system('gunzip '+args.clean_fastq_1+";"+'gunzip '+args.clean_fastq_2)
a = os.popen("cat "+clean_fastq_1 +"| wc -l").read()
fastqreads=int(int(a)/2)
print(fastqreads)
if os.path.exists(clean_fastq_1):
    print('*********BWA starting :', time.asctime(time.localtime(time.time())),'\n')
    os.system('bwa index -a bwtsw -p '+virusfasta+' ' +virusfasta+ "> /dev/null 2>&1")
    os.system('bwa mem -t '+threads+' '+virusfasta+' ' +clean_fastq_1+' '+clean_fastq_2+'| samtools sort |samtools view -bh -f 2 -o '+virus_sorted_bam+ " > /dev/null 2>&1")
    os.system('samtools mpileup -f '+virusfasta+' '+virus_sorted_bam+'  > '+pipeup)
    os.system('python '+softdir+'/count.py '+pipeup+' '+consensus)
    os.system('rm '+pipeup)
    print('*********MetaCompass starting :', time.asctime(time.localtime(time.time())))
    os.system('python3  '+MetaCompass_dir+'/go_metacompass.py -r '+virusfasta+' -1 '+clean_fastq_1+' -2 '+clean_fastq_2+' -l '+readslen+' -o '+metacompassout+' -t '+threads)
    if os.path.exists(genome):
        os.system("sed '/len=/,+1d' "+genome+' > '+genome1)
    elif os.path.exists(contig):
        os.system("sed '/len=/,+1d' "+contig+' > '+contig1)
    print('*********RagTag starting :', time.asctime(time.localtime(time.time())))
    os.system("mkdir -p "+ragtag_output)
    metacompass_genomes=open(outdir+'/metacompass/metacompass_output/metacompass.genomes_coverage.txt','r')
    data = pd.read_csv(metacompass_genomes,header=0,sep='\t')
    all=[]
    all=data["ref_ID"].tolist()
    os.system("python "+softdir+"/seperate.py "+virus_fasta+" "+outdir)
    if os.path.exists(genome):
        for id in all:
            os.system("grep "+id+" -A 1 "+genome1+" > "+metacompassout+'/'+id+'.metacompass.fa')
            os.system("ragtag.py scaffold  "+virusfasta+' '+metacompassout+'/'+id+'.metacompass.fa  --debug -w -o '+ragtag_output+'/'+id+'  -t 20 -u -C')
            os.system("sed '/Chr0/,+1d' "+ragtag_output+"/"+id+"/ragtag.scaffold.fasta >  "+ragtag_output+"/"+id+"/ragtag.scaffold.fasta1")
            os.system("echo virus: "+id+" Quast the Genome fraction:")
            os.system("metaquast.py "+ragtag_output+"/"+id+"/ragtag.scaffold.fasta1  -r  "+outdir+'/ref/'+id+'.fa'+" --output-dir "+outdir+"/metaquast/"+id+" 1>"+outdir+"/metaquast.log")
            os.system("python  "+softdir+"/result.py "+id+' '+outdir+"/metaquast/"+id+"/summary/TXT/Genome_fraction.txt "+str(fastqreads)+ ' '+outdir)
            os.system("cat "+outdir+"/result/"+id+"/result.txt")
    elif os.path.exists(contig):
        for id in all:
            os.system("grep "+id+" -A 1 "+contig1+" > "+metacompassout+'/'+id+'.metacompass.fa')
            os.system("ragtag.py scaffold  "+virusfasta+' '+metacompassout+'/'+id+'.metacompass.fa')
            os.system("ragtag.py scaffold  "+virusfasta+' '+metacompassout+'/'+id+'.metacompass.fa  --debug -w -o '+ragtag_output+'/line  -t 20 -u -C')
            os.system("sed '/Chr0/,+1d' "+ragtag_output+"/"+id+"/ragtag.scaffold.fasta >  "+ragtag_output+"/"+id+"/ragtag.scaffold.fasta1")

os.system('samtools depth '+virus_sorted_bam+' > '+outdir+'/result/depth.txt')
os.system("python  "+softdir+"/pic.py "+outdir+'/result/depth.txt '+outdir+'/result/')
os.system('rm '+outdir+'/result/depth.txt')
print('Finish time:', time.asctime(time.localtime(time.time())))