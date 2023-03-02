# VIGA

VIGA is an effective workflow for Virus Identification and Genome Assembly from NGS data

## Installation

Download VIGA with Git from GitHub

```
git clone https://github.com/viralInformatics/VIGA.git
```

### Software

```
conda env create -n VIGA -f VIGA.yml
```

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
bash $file/runall.sh -i sample -t 30
```

#### 2 Step by step

cd test_data 

put paired-end gz in ./test_data/fastq/

##### Step1: Virus Identification

Input fastq file must be paired-end sequencing, put it under the ./test_data/fastq/ file, and the sample name is the name before the underline of _1.fastq.gz

```
Usage: bash $file/bin/0_run1.sh [-h] [-i input sample name ] [-t threads]
      -h  <Display this help message.>
      -i   <input sample name. eg ERR3276898_1.fastq.gz sample name is ERR3276898>
      -t   <threads number. eg 30>

eg. bash $file/bin/0_run1.sh -i ERR3276898 -t 30


Required parameter:

    -i    < input sample name >
    -t    < threads number >
```

##### Step2: Virus Genome Assemble

The result of the first step of identification is located in speciesfinal.txt under the ./test_data/classify/ folder, and the sequence file is located in ./test_data/genome/virus/sample.fa

```
Usage: bash $file/bin/0_run2.sh [-h] [-i input sample name ] [-t threads]
      -h  <Display this help message.>
      -i   <input sample name. eg ERR3276898_1.fastq.gz sample name is ERR3276898>
      -t   <threads number. eg 30>

eg. bash $file/bin/0_run2.sh -i ERR3276898 -t 30


Required parameter:

    -i    < input sample name >
    -t    < threads number >
```

The final result is under the ./test_data/genome/sample/result folder
