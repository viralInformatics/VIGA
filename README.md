![image-20230316212835572](img/README/image-20230316212835572.png)

# VIGA

VIGA is an effective workflow for Virus Identification and Genome Assembly from NGS data

## Installation

### Step1: Download VIGA

Download VIGA with Git from GitHub

```
git clone https://github.com/viralInformatics/VIGA.git
```

or Download ZIP to local

### Step 2: Download Database

```
1. download taxdmp.zip [Index of /pub/taxonomy (nih.gov)](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/) and unzip taxdmp.zip and put it in ./db/

2. download "prot.accession2taxid" file from https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/

3. download "RefSeqVirusProtein" file from
wget -c ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.protein.faa.gz
gzip -d viral.1.protein.faa.gz
mv viral.1.protein.faa RefSeqVirusProtein

4. download "nr" file from
wget -c ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz
or ascp -T  -i  asperaweb_id_dsa.openssh --host=ftp.ncbi.nih.gov --user=anonftp --mode=recv /blast/db/FASTA/nr.gz ./
gzip -d nr.gz

5. Use Diamond v2.0.11.149 to create two separate databases as the indexing libraries in the current version are incompatible with each other.

6. In order to set up a reference database for DIAMOND, the makedb command needs to be executed with the following command line:
diamond makedb --in YourPath/RefSeqVirusProtein  -d Diamond_RefSeqVirusProtein --taxonmap YourPath/prot.accession2taxid --taxonnodes YourPath/nodes.dmp
diamond makedb --in nr -d Dimond_nr --taxonmap YourPath/prot.accession2taxid --taxonnodes YourPath/nodes.dmp

```

### Step 3: Installation of dependent software

#### Installing Some Software Using Conda

```
conda install fastp=0.12.4 trinity=2.8.5 diamond=2.0.11.149 ragtag=2.1.0 quast=5.0.2
```

#### Manual Installation of MetaCompass

https://github.com/marbl/MetaCompass

### Step 4: Python Dependencies

Base on python 3.6.8

```
pip install pandas=1.1.5 numpy=1.19.5  matplotlib=3.3.4  biopython=1.79
```



## Usage

Note: If the reference genomes of viruses are known, please place the sequences in the ./Ref folder and name them after the sample using the format SRAID.fa, such as ERR3253399.fa, then proceed to run **Step 2: Virus Genome Assembly** directly.

### Step1: Virus Identification

```
#Paire-ended files:
usage: 0_run1_paired.py [-h] [--evalue EVALUE] --fastq_1 FASTQ_1 --fastq_2
                        FASTQ_2 --outdir OUTDIR --Diamond_VirusProtein_db
                        DIAMOND_VIRUSPROTEIN_DB --Diamond_nr_db DIAMOND_NR_DB
                        [--threads THREADS]
eg. 
python /data/12T/fp/software/VIGA/bin/0_run1_paired.py --fastq_1 /data/12T/fp/plant_1000/testVIGA/fastq/ERR3253399_1.fastq.gz  --fastq_2  /data/12T/fp/plant_1000/testVIGA/fastq/ERR3253399_2.fastq.gz --outdir  /data/12T/fp/plant_1000/testVIGA/test --Diamond_VirusProtein_db /data/12T/fp/plant_1000/test_mock/db/diamondRefSeqVirusProtein --Diamond_nr_db /data/12T/fp/plant_1000/test_mock/db/diamond-nr --evalue 1e-5  --threads 10 

#Single-end sequencing files only support step1 virus identification
usage: 0_run1_single.py [-h] [--evalue EVALUE] --fastq FASTQ --outdir OUTDIR --Diamond_VirusProtein_db DIAMOND_VIRUSPROTEIN_DB --Diamond_nr_db DIAMOND_NR_DB [--threads THREADS]

eg. 
python SoftwarePath/0_run1_single.py --fastq YourPath/sample.fastq.gz --outdir  YourPath/test --Diamond_VirusProtein_db YourPath/db/diamondRefSeqVirusProtein --Diamond_nr_db YourPath/db/diamond-nr --evalue 1e-5  --threads 10 
```

### Step2: Virus Genome Assemble

The result of the first step of identification is located in speciesfinal.txt under the ./test/Classify/ folder, and the sequence file is located in ./test/Ref/sample.fa



```
usage: 0_run2.py [-h] [--len LEN] --clean_fastq_1 CLEAN_FASTQ_1
                 --clean_fastq_2 CLEAN_FASTQ_2 --virus_fasta VIRUS_FASTA
                 --outdir OUTDIR --MetaCompass_dir METACOMPASS_DIR
                 [--threads THREADS]
0_run2.py: the following arguments are required: --clean_fastq_1, --clean_fastq_2, --virus_fasta, --outdir, --MetaCompass_dir

eg. 
python /data/12T/fp/software/VIGA/bin/0_run2.py --clean_fastq_1 /data/12T/fp/plant_1000/testVIGA/test/Fastp/ERR3253399_1.clean.fastq.gz  --clean_fastq_2 /data/12T/fp/plant_1000/testVIGA/test/Fastp/ERR3253399_2.clean.fastq.gz --virus_fasta /data/12T/fp/plant_1000/testVIGA/test/Ref/ERR3253399.fa --outdir /data/12T/fp/plant_1000/testVIGA/test/Genome --MetaCompass_dir /data/12T/fp/software/MetaCompass --threads 20
```

The final result is under the ./test/Genome/result folder, the depth picture is in result/picture folder

| *Column*    | *Description*                                    |
| ----------- | ------------------------------------------------ |
| VIGA(%)     | Viral Genome Completeness Evaluated by MetaQuast |
| abundance   | FPKM                                             |
| coverage(%) | Broadness of the virus genome                    |
| depth_cov   | Depth of the covered the viral genome            |
| depth_all   | Depth of all the virus genome                    |
