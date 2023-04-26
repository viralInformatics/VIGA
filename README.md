![image-20230316212835572](img/README/image-20230316212835572.png)

# VIGA

VIGA is an effective workflow for Virus Identification and Genome Assembly from NGS data

## Installation

Download VIGA with Git from GitHub

```
git clone https://github.com/viralInformatics/VIGA.git
```

### Database

```
1. download taxdmp.zip [Index of /pub/taxonomy (nih.gov)](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/) and unzip taxdmp.zip and put it in ./db/
2. download "prot.accession2taxid" file from https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/
3. Use Diamond v2.0.11.149 to create two separate databases as the indexing libraries in the current version are incompatible with each other.

diamond makedb --in YourPath/RefSeqVirusProtein  -d Diamond_RefSeqVirusProtein --taxonmap YourPath/prot.accession2taxid --taxonnodes YourPath/nodes.dmp
diamond makedb --in nr -d Dimond_nr --taxonmap YourPath/prot.accession2taxid --taxonnodes YourPath/nodes.dmp
```

### Software

```
conda env create -n VIGA -f VIGA.yml
```

#### MetaCompass

https://github.com/marbl/MetaCompass

### Python Dependencies

Base on python 3.6.8

```
pip install -r requirements.txt
```

## Usage

### Running VIGA pipeline

#### 1 One step

If the library of diamon blastx has been built in the corresponding db folder, you can run it directly

```
usage: 0_runall.py [-h] [--evalue EVALUE] --fastq_1 FASTQ_1 --fastq_2 FASTQ_2
                   --outdir OUTDIR --Diamond_VirusProtein_db
                   DIAMOND_VIRUSPROTEIN_DB --Diamond_nr_db DIAMOND_NR_DB
                   [--threads THREADS] --virus_fasta VIRUS_FASTA
                   --clean_fastq_1 CLEAN_FASTQ_1 --clean_fastq_2 CLEAN_FASTQ_2
                   --MetaCompass_dir METACOMPASS_DIR
0_runall.py: the following arguments are required: --fastq_1, --fastq_2, --outdir, --Diamond_VirusProtein_db, --Diamond_nr_db, --virus_fasta, --clean_fastq_1, --clean_fastq_2, --MetaCompass_dir

eg. python SoftwarePath/0_runall.py  --clean_fastq_1 YourPath/test/Fastp/ERR3253399_1.clean.fastq.gz  --clean_fastq_2 YourPath/test/Fastp/ERR3253399_2.clean.fastq.gz   --virus_fasta YourPath/test/Ref/ERR3253399.fa --outdir YourPath/test/ --MetaCompass_dir YourPath/software/MetaCompass --threads 20   --fastq_1 YourPath/ERR3253399_1.fastq.gz  --fastq_2  YourPath/ERR3253399_2.fastq.gz  --Diamond_VirusProtein_db YourPath/db/diamondRefSeqVirusProtein --Diamond_nr_db YourPath/db/diamond-nr --evalue 1e-5
```

#### 2 Step by step

##### Step1: Virus Identification

```
usage: 0_run1_paired.py [-h] [--evalue EVALUE] --fastq_1 FASTQ_1 --fastq_2
                        FASTQ_2 --outdir OUTDIR --Diamond_VirusProtein_db
                        DIAMOND_VIRUSPROTEIN_DB --Diamond_nr_db DIAMOND_NR_DB
                        [--threads THREADS]
0_run1_paired.py: error: argument --fastq is required

eg. 
python SoftwarePath/0_run1_paired.py --fastq_1 YourPath/ERR3253399_1.fastq.gz  --fastq_2  YourPath/ERR3253399_2.fastq.gz --outdir  YourPath/test --Diamond_VirusProtein_db YourPath/db/diamondRefSeqVirusProtein --Diamond_nr_db YourPath/db/diamond-nr --evalue 1e-5  --threads 10 


#Single-ended files can only be identified for Step1

usage: 0_run1_single.py [-h] [--evalue EVALUE] --fastq FASTQ --outdir OUTDIR --Diamond_VirusProtein_db
                        DIAMOND_VIRUSPROTEIN_DB --Diamond_nr_db DIAMOND_NR_DB
                        [--threads THREADS]
0_run1_single.py: error: argument --fastq is required

eg. 
python SoftwarePath/0_run1_paired.py --fastq_1 YourPath/ERR3253399_1.fastq.gz  --fastq_2  YourPath/ERR3253399_2.fastq.gz --outdir  YourPath/test --Diamond_VirusProtein_db YourPath/db/diamondRefSeqVirusProtein --Diamond_nr_db YourPath/db/diamond-nr --evalue 1e-5  --threads 10 
```

##### Step2: Virus Genome Assemble

The result of the first step of identification is located in speciesfinal.txt under the ./test/Classify/ folder, and the sequence file is located in ./test/Ref/sample.fa

```
usage: 0_run2.py [-h] [--len LEN] --clean_fastq_1 CLEAN_FASTQ_1
                 --clean_fastq_2 CLEAN_FASTQ_2 --virus_fasta VIRUS_FASTA
                 --outdir OUTDIR --MetaCompass_dir METACOMPASS_DIR
                 [--threads THREADS]
0_run2.py: the following arguments are required: --clean_fastq_1, --clean_fastq_2, --virus_fasta, --outdir, --MetaCompass_dir

eg. 
python SoftwarePath/0_run2.py --clean_fastq_1 YourPath/Fastp/ERR3253399_1.clean.fastq.gz  --clean_fastq_2 YourPath/Fastp/ERR3253399_2.clean.fastq.gz --virus_fastaYourPath/test/Ref/ERR3253399.fa --outdir /YourPath/test/Genome --MetaCompass_dir YourPath/software/MetaCompass --threads 20
```

The final result is under the ./test/Genome/result folder; the depth picture is in result/picture folder
