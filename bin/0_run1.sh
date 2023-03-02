#!/usr/bin/env bash
#source /home/root640/software/miniconda3/bin/activate fp
#
#name=$1
file=$(pwd)
echo -e "##############################################################\n"
#echo -e $(date +%F%n%T)'------Script start: \n'
echo  -e 'The working directory is : '$file'\n'
echo -e "##############################################################\n"

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

name=""
nThread=""
# 参数传递


while getopts ":hit:" opt; do
    case ${opt} in
        h ) print_help;;
        i ) name="$OPTARG";;
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


echo "Processing file: $name"

if [ $# = 0 ]
then
    helpdoc
    exit 1
fi

# file
mkdir -p $file/fastq;mkdir -p $file/trinity;mkdir -p $file/blastvirus;mkdir -p $file/blastvirus_fa;mkdir -p $file/blastnr;mkdir -p $file/classify;mkdir -p $file/genome/virus;mkdir -p $file/log

path=$0
software=$(dirname $(dirname $path))

################################
########arguments############
################################
######################################################
for sra in $name;do
sra1=${sra##*/}
name=${sra1%%.*}
#fq_1_gz=$file/fastq/$name'_1.fastq.gz'
fq_1_gz=`ls $file/fastq/$name*1*gz|head -1`
fq_2_gz=`ls $file/fastq/$name*2*gz|head -1`
#fq_2_gz=$file/fastq/$name'_2.fastq.gz'
fq_gz=$file/fastq/$name'.fastq.gz'
clean_fq_1_gz=$file/fastq/$name'_1.clean.fastq.gz' 
clean_fq_2_gz=$file/fastq/$name'_2.clean.fastq.gz' 
clean_fq_gz=$file/fastq/$name'.clean.fastq.gz'
clean_fq_1=$file/fastq/$name'_1.clean.fastq' 
clean_fq_2=$file/fastq/$name'_2.clean.fastq' 
clean_fq=$file/fastq/$name'.clean.fastq'

time1=$(date +%F%n%T)
date1=$(date +"%Y-%m-%d %H:%M:%S")
sys_date1=$(date -d "$date1" +%s)
########core###################
echo -e "##############################################################\n"
echo -e $(date +%F%n%T)'------fastp start: \n'
echo -e "##############################################################\n"
echo -e "fastp command: fastp -i $fq_1_gz  -I  $fq_2_gz  -o  $clean_fq_1   -O $clean_fq_2  -j $file/fastq/$name.json -h $file/fastq/$name.html -w $nThread'\n'"
if [  ! -f $clean_fq_1  ];then
   fastp -i $fq_1_gz  -I  $fq_2_gz  -o  $clean_fq_1   -O $clean_fq_2  -j $file/fastq/$name.json -h $file/fastq/$name.html -w $nThread
elif [  -f $fq_1_gz   ];then
   echo 'fastp have done!'
fi

sed -i 's/ //g' $clean_fq_1
sed -i 's/ //g' $clean_fq_2

echo -e "##############################################################\n"
echo -e $(date +%F%n%T)'------Trinity start: \n'
echo -e "##############################################################\n"

if [  -f $clean_fq_1  ];then
	echo 'Trinity command: Trinity --seqType fq --left  '$clean_fq_1' --right '$clean_fq_2' --output '$file'/trinity/trinity_out_dir_'$name' --CPU 26 --max_memory 50G';Trinity --seqType fq --left $clean_fq_1 --right $clean_fq_2 --output $file/trinity/trinity_out_dir_$name --CPU 26 --max_memory 50G
elif [ ! -f $clean_fq_1  ];then
	echo -e 'Clean fastq file not exists under ./fastq/ dir!.'
fi

if [ ! -f  $file/trinity/trinity_out_dir_$name/Trinity.fasta ];then
    echo 'Trinity failed, Megahit starting......';megahit -1  $clean_fq_1   -2  $clean_fq_2 -o $file/trinity/trinity_out_dir_$name/megahit;mv $file/trinity/trinity_out_dir_$name/megahit/final.contigs.fa $file/trinity/trinity_out_dir_$name/Trinity.fasta
fi

echo -e "##############################################################\n"
echo -e $(date +%F%n%T)'------Trinity end: \n'
echo -e "##############################################################\n"

echo -e "##############################################################\n"
echo -e $(date +%F%n%T)'------Diamond blastx virusrefseq start : \n'
echo -e "##############################################################\n"

diamond blastx -q $file/trinity/trinity_out_dir_$name/Trinity.fasta  --db /data/12T/fp/plant_1000/test_mock/db/diamondRefSeqVirusProtein -o $file/blastvirus/$name'.vp.txt'  -e 1e-5   --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore  slen  stitle salltitles qcovhsp nident staxids 1>$file/log/diamond.log 2>$file/log/diamond.error

echo -e "##############################################################\n"
echo -e $(date +%F%n%T)'------Get fasta start : \n'
echo -e "##############################################################\n"

cat $file/blastvirus/$name'.vp.txt'  |awk '{print $1}'|uniq > $file/blastvirus/$name'vp1nd'

cat $file/trinity/trinity_out_dir_$name/Trinity.fasta | seqkit grep -f  $file/blastvirus/$name'vp1nd' > $file/blastvirus_fa/$name'.fa'

echo -e "##############################################################\n"
echo -e $(date +%F%n%T)'------Diamond blastx nr start : \n'
echo -e "##############################################################\n"

diamond blastx -q $file/blastvirus_fa/$name'.fa'  --db /data/12T/fp/plant_1000/test_mock/db/diamond-nr -o $file/blastnr/$name'_diamond.txt'  -e 1e-5   --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore  slen  stitle salltitles qcovhsp nident staxids 1>$file/log/diamond.log 2>$file/log/diamond.error

cat $file/blastnr/$name'_diamond.txt' | sort -k1,1 -k12,12gr -k11,11g -k3,3gr | sort -u -k1,1 --merge> $file/blastnr/$name'_diamond.besthit'

python $software/bin/filter.py $file/blastnr/$name'_diamond.besthit'  $file/blastnr/$name'_virussure.besthit'    $file/classify/$name'_species' $file/classify/$name'_genus' $file/classify/$name'_family' $name $software

rm  $file/blastvirus/$name'vp1nd'

time2=$(date +%F%n%T)
date2=$(date +"%Y-%m-%d %H:%M:%S")
sys_date2=$(date -d "$date2" +%s)
time=`expr $sys_date2 - $sys_date1`
echo -e "##############################################################\n"
echo -e "The run1.sh total waste time : $time"s"\n"
echo -e "##############################################################\n"


done
