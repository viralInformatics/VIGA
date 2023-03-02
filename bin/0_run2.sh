#!/bin/bash
#author fuping
#2022/3/7
##source /home/root640/software/miniconda3/bin/activate fp
file=$(pwd)
time1=$(date +%F%n%T)
echo $time1 'run2 start......'
date1=$(date +"%Y-%m-%d %H:%M:%S")
sys_date1=$(date -d "$date1" +%s)

virusid=$1
name=$1
################################
########help###################
########help###################

function print_help {
    cat <<EOF
Description:

   This is a help document
 * This script is used for virus identification *

Usage:

Usage: $0 [-h]
      -h  <Display this help message.>

eg. bash $0 -i ERR3276898 -t 30

      -i   <input sample name. eg ERR3276898_1.fastq.gz sample name is ERR3276898>  
      -t   <threads number. eg 30>

Required parameter:

    -i    < input sample name >
    -t    < threads number >

EOF
}

SRAfile=""
nThread=""
# 参数传递


while getopts ":hit:" opt; do
    case ${opt} in
        h ) print_help;;
        i )  name="$OPTARG";;
        t ) nThread="$OPTARG";;
        \? ) echo "Invalid option: -$OPTARG" 1>&2; exit 1;;
        : ) echo "Option -$OPTARG requires an argument." 1>&2; exit 1;;

    esac
done

if [ $# = 0 ]
then
    print_help
    exit 1
fi

if [ -z "$name" ]; then
  echo "Error: No input sample specified." >&2
  exit 1
fi
if [ -z "$nThread" ]; then
  echo "Error: No nThread specified." >&2
  exit 1
fi

echo -e "##############################################################\n"
echo -e  "sample: "$name"\n"
echo -e "##############################################################\n"
################################
########arguments############
################################

mkdir -p $file/fastq;mkdir -p $file/trinity;mkdir -p $file/blastvirus;mkdir -p $file/blastvirus_fa;mkdir -p $file/blastnr;mkdir -p $file/classify;mkdir -p $file/genome/virus;mkdir -p $file/log

if [ -f $file/genome/$name ];then
	rm -rf $file/genome/$name
elif [ ! -f $file/genome/$name ];then
	mkdir -p $file/genome/$name
fi
virusid_fa=$file/genome/virus/$virusid'.fa'
echo -e 'Reference virus in the '$virusid_fa'\n'
clean_fq_1_gz=$file/fastq/$name'_1.clean.fastq.gz' 
clean_fq_2_gz=$file/fastq/$name'_2.clean.fastq.gz' 
clean_fq_gz=$file/fastq/$name'.clean.fastq.gz'
clean_fq_1=$file/fastq/$name'_1.clean.fastq' 
clean_fq_2=$file/fastq/$name'_2.clean.fastq' 
clean_fq=$file/fastq/$name'.clean.fastq'
virusid_sorted_bam=$file/genome/$name/$virusid'.sorted.bam'

path=$0
software=$(dirname $(dirname $path))

if [ -f $clean_fq_1_gz ];then
	echo 'gunzip '$clean_fq_1_gz' and '$clean_fq_2_gz;gunzip $clean_fq_1_gz;gunzip $clean_fq_2_gz
elif [ -f $clean_fq_1 ];then
	echo -e 'Check '$clean_fq_1';'$$clean_fq_2' done!\n'
fi

################################
################################
echo -e "##############################################################\n"
echo -e $(date +%F%n%T)'------bwa mem start: \n'
echo -e "##############################################################\n"
#####
bwa index -a bwtsw -p $virusid_fa $virusid_fa

if [  -f $virusid_sorted_bam ];then
    size=`wc -c $virusid_sorted_bam | awk '{print $1}'`
    if [ $size -le 100 ];then
        rm -f $virusid_sorted_bam
    fi
fi

if [ ! -f $virusid_sorted_bam ];then
    bwa mem -t 20 $virusid_fa  $clean_fq_1 $clean_fq_2  | samtools sort |samtools view -bh -f 2 -o $virusid_sorted_bam 1>$file/log/bwa.log 2>$file/log/bwa.error;samtools index $virusid_sorted_bam
fi

if [ -f $virusid_sorted_bam ];then
    samtools mpileup -f $file/genome/virus/$name'.fa'  $file/genome/$name/$file'.sorted.bam' > $file/genome/$name/out.pipeup;python count.py $file/genome/$name/out.pipeup $file/genome/$name/consensus.fa;rm $file/genome/$name/out.pipeup
fi
#####MetaCompass
echo -e "##############################################################\n"
echo -e $(date +%F%n%T)'-------MetaCompass  start: \n'
echo -e "##############################################################\n"


#if  [  -d $file/genome/$name/metacompass ];then
 #   rm -rf $file/genome/$name/metacompass
#fi

#if users install by themselves, the path is :./MetaCompass/go_metacompass.py

echo python3 /data/12T/fp/software/MetaCompass/go_metacompass.py -r $virusid_fa  -1 $clean_fq_1 -2 $clean_fq_2 -l 100 -o $file/genome/$name/metacompass -t 30
python3 /data/12T/fp/software/MetaCompass/go_metacompass.py -r $virusid_fa  -1 $clean_fq_1 -2 $clean_fq_2 -l 100 -o $file/genome/$name/metacompass -t 30 

genome=$file/genome/$name/metacompass/metacompass_output/metacompass.final.ctg.fa
contig=$file/genome/$name/assembly/contigs.fasta

genome1=$file/genome/$name/metacompass/metacompass_output/metacompass.final.ctg.fa1
contig1=$file/genome/$name/assembly/contigs.fasta1


if  [  -f  $genome  ];then
    sed '/len=/,+1d'  $file/genome/$name/metacompass/metacompass_output/metacompass.final.ctg.fa > $genome1
elif  [ -f  $contig  ];then
    sed '/len=/,+1d'  $file/genome/$name/metacompass/assembly/contigs.fasta  > $contig1
fi

echo -e "##############################################################\n"
echo -e $(date +%F%n%T)'-------ragtag  start: \n'
echo -e "##############################################################\n"

mkdir -p $file/genome/$name/ragtag_output
chr=$(cat $file/genome/$name/metacompass/metacompass_output/metacompass.genomes_coverage.txt |awk '{ if (NR >1 ){print $1}}'|tr '\n' ',')

for id in ${chr//,/ };do
    if  [  -f  $genome  ];then
        grep $id -A 1 $genome1 > $file/genome/$name/metacompass/metacompass_output/$id'.metacompass.fa';ragtag.py scaffold  $file/genome/virus/$id'.fa'  $file/genome/$name/metacompass/metacompass_output/$id'.metacompass.fa'  --debug -w -o $file/genome/$name/ragtag_output/$id  -t 20 -u -C;sed '/Chr0/,+1d' $file/genome/$name/ragtag_output/$id/ragtag.scaffold.fasta >  $file/genome/$name/ragtag_output/$id/ragtag.scaffold.fasta1
    elif  [ -f  $contig  ];then
         grep $id -A 1 $contig1 >  $contig1'.fa';ragtag.py scaffold  $file/genome/virus/$id'.fa'  $contig1'.fa' --debug -w  -o $file/genome/$name/ragtag_output/$id  -t 20 -u -C;sed '/Chr0/,+1d' $file/genome/$name/ragtag_output/$id/ragtag.scaffold.fasta >  $file/genome/$name/ragtag_output/$id/ragtag.scaffold.fasta1
    fi
done


echo  -e  $(date +%F%n%T)'------get trinity fasta......\n'
trinity_fa=$file/blastvirus_fa/$name'.fa'
echo -e "##############################################################\n"
echo -e  $(date +%F%n%T)'------ragtag plus trinity fasta start......\n'
echo -e "##############################################################\n"

allfa=$file/genome/$name/all.fa
if [  -f $allfa  ];then
   rm $allfa
elif [ ! -f $allfa ];then
   echo $allfa'not exists! Get at '$(date +%F%n%T)
fi

cat $trinity_fa >> $allfa

if  [  -f  $genome  ];then
    cat $file/genome/$name/metacompass/metacompass_output/metacompass.final.ctg.fa1 >>  $allfa
elif  [ -f  $contig  ];then
    cat $file/genome/$name/metacompass/assembly/contigs.fasta1  >>  $allfa
fi

cat  $file/genome/$name/ragtag_output/*/ragtag.scaffold.fasta1 >>  $allfa
ragtag.py scaffold  $virusid_fa  $allfa  --debug -w  -o $file/genome/$name/ragtag_output_trinity  -t 20 -u -C

sed '/Chr0/,+1d' $file/genome/$name/ragtag_output_trinity/ragtag.scaffold.fasta >  $file/genome/$name/ragtag_output_trinity/ragtag.scaffold.fasta1

echo -e "##############################################################\n"
echo -e $(date +%F%n%T)'------metaquast is start......\n'
echo -e "##############################################################\n"

b=2
if [  -f $clean_fq_1_gz ];then
    a=$(zcat $clean_fq_1_gz|wc -l);fastqreads=`expr $a / $b`
elif [ -f $clean_fq_1 ];then
    a=$(cat $clean_fq_1|wc -l);fastqreads=`expr $a / $b`   #fastqreads=`expr $(cat $clean_fq_1|wc -l)/$chushu`
fi
echo 'FastqReads number: '$fastqreads




for id in ${chr//,/ };do
    if [ ! -f $file/genome/virus/$id'.fa' ];then
        python $software/bin/seperate.py $name
    fi
	echo 'virus: '$id
    metaquast.py $file/genome/$name/metacompass/metacompass_output/$id'.metacompass.fa'  -r $file/genome/virus/$id'.fa'  --output-dir $file/genome/$name/metaquast/$id/metacompass 1>$file/log/metaquast.log 2>$file/log/metaquast.error
    metaquast.py $file/genome/$name/ragtag_output/$id/ragtag.scaffold.fasta1 -r   $file/genome/virus/$id'.fa'  --output-dir $file/genome/$name/metaquast/$id/ragtag_output 1>$file/log/metaquast.log 2>$file/log/metaquast.error
    metaquast.py  $file/genome/$name/ragtag_output_trinity/ragtag.scaffold.fasta   -r  $file/genome/virus/$id'.fa' --output-dir $file/genome/$name/metaquast/$id/ragtag_output_trinity 1>$file/log/metaquast.log 2>$file/log/metaquast.error
    metaquast.py  $file/trinity/'trinity_out_dir_'$name/Trinity.fasta   -r  $file/genome/virus/$id'.fa'  --output-dir $file/genome/$name/metaquast/$id/trinity 1>$file/log/metaquast.log 2>$file/log/metaquast.error
    meta=$file/genome/$name/metaquast/$id/metacompass/summary/TXT/Genome_fraction.txt
    rag_fra=$file/genome/$name/metaquast/$id/ragtag_output/summary/TXT/Genome_fraction.txt
    rag_trinity_fra=$file/genome/$name/metaquast/$id/ragtag_output_trinity/summary/TXT/Genome_fraction.txt
    trinity=$file/genome/$name/metaquast/$id/trinity/summary/TXT/Genome_fraction.txt
    echo -e "############################  $id details :#############################\n"
    python $software/bin/result.py $meta $rag_fra $rag_trinity_fra $trinity $fastqreads $id $name
    echo -e "############################ $id result :#############################\n"
    cat $file/genome/$name/result/$id/result.txt
    echo -e  "\n######################################################################\n"
done

echo -e "##############################################################\n"
echo -e $(date +%F%n%T)'------new virus dicovery is start......\n'
echo -e "##############################################################\n"
python $software/bin/getLCA.py $name $software

time2=$(date +%F%n%T)

date2=$(date +"%Y-%m-%d %H:%M:%S")
sys_date2=$(date -d "$date2" +%s)
time=`expr $sys_date2 - $sys_date1`
echo -e "##############################################################\n"
echo -e "run2.sh total waste time : $time"s"\n"
echo -e "##############################################################\n"

