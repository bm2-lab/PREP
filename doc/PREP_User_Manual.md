# PREP User Manual 


## Table of Contents
1. [General Description](#general-description)  
2. [Dependencies](#dependencies)  
    - [Required software](#required-software)  
    - [Python packages](#python-packages) 
3. [Installation via Docker](#installation-via-docker)  
4. [Installation from source](#installation-from-source)  
5. [Usage](#usage)  
6. [Input Files](#input-files)  
    - [Input Files](#input-files (required))  
7. [Setting parameters](#seting-parameters) 
8. [Output Files](#output-files)  
    - [Column explanation](#column-explanation)  
9. [Contact](#contact)
10. [Algorithmic Flow Chart](#algorithmic-flow-chart)

## General Description

RNA editing is a source of transcriptomic diversity, mainly in non-coding regions, and is altered in cancer. Recent studies demonstrated that A-to-I RNA editing events are manifested at the proteomic level and contribute to protein heterogeneity in cancer. Given somatic RNA-editing mutation as input, PREP identify and evaluate the potential immunogenicity of RNA editing based peptides. Detailed information please refer to citation.

## Dependencies  

#### Hardware:
PREP currently tested on x86_64 on ubuntu 16.04.

#### Required software:
* [Python 2.7](https://www.python.org/downloads/release/python-2712/)
* [NetMHCpan 4.0](http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?netMHCpan)
* [Variant Effect Predictor (VEP)](https://github.com/Ensembl/ensembl-vep)
* [BWA](https://github.com/lh3/bwa)
* [STAR](https://github.com/alexdobin/STAR)
* [samtools](https://github.com/samtools)
* [Optitype](https://github.com/FRED-2/OptiType)
* [GATK 3.8](https://software.broadinstitute.org/gatk/best-practices/)
* [Picard tools](https://broadinstitute.github.io/picard/)
* [Java 8](https://java.com/en/download/help/linux_x64rpm_install.xml)
* [kallisto](http://pachterlab.github.io/kallisto/)
* [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
* [vcftools](http://vcftools.sourceforge.net/)
* [blast](http://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
* [tabix](http://www.htslib.org/doc/tabix.html)
* [gawk]()

#### Required Python package:
* [yaml](https://pypi.org/project/yaml-1.3/)
* [XGboost](https://pypi.org/project/xgboost/)
* [biopython](https://pypi.org/project/biopython/)
* [scikit-learn==0.19.1](https://pypi.org/project/scikit-learn/)
* [pandas](https://pypi.org/project/pandas/)
* [numpy](https://pypi.org/project/numpy/)
* [Pyomo](https://pypi.org/project/Pyomo/)
* [tables](https://pypi.org/project/tables/)
* [pysam](https://pypi.org/project/pysam/)
* [future](https://pypi.org/project/future/)
* [multiprocessing](https://pypi.org/project/multiprocessing/)
* [subprocess](https://pypi.org/project/subprocess/)
* [math](https://pypi.org/project/math/)
* [matplotlib](https://pypi.org/project/matplotlib/)


## Installation via Docker
Docker image of PREP is at https://hub.docker.com/r/bm2lab/prep/.

1. Install Docker on your computer and make sure it works.

2. Call docker `pull bm2lab/prep` which will download the Docker image.

3. Run the image in interactive mode with your dataset:
        
		docker run -it -v /your/path/to/dataset/:/home/bioworker/dataset bm2lab/prep /bin/bash

4. Change directory into /home/bioworker/project/PREP:

		cd /home/bioworker/project/PREP

5. Download reference data:

		bash data_download.sh

6. Edit `config.yaml` and fill the proper path of input files.

7. Run the program with follow commands:

		python PREP.py RE -i config.yaml


## Installation from source

1. Install all software, python packages and R packages listed above, and make sure each software and package works in your system. 
2. Install multiprocessing and other packages with the `pip` command:

        pip install -U multiprocessing
        ...
 
4. Download or clone the PREP repository to your local system:

        git clone https://github.com/bm2-lab/PREP.git

5. Reference data includes genome fasta, peptide(GRCh38 build) could be downloaded and processed through:

        bash data_download.sh
        
    a few reference data would be in the fold `database` and processed by custom script in order to run the pipeline, including:

        [Fasta] 
	
        This fold contains the reference fasta file, its bwa index and some other files result from `huamn.fasta`:
        Homo_sapiens.GRCh38.dna.primary_assembly.fa
        Homo_sapiens.GRCh38.dna.primary_assembly.fa.amb	
        Homo_sapiens.GRCh38.dna.primary_assembly.fa.ann	
        etc...
	
        [Annotation file] 
	
        This fold contains the vcf file used to run RNAEditor:
        1000GenomeProject.vcf
        HAPMAP.vcf
        ESP.vcf
        dbSNP.vcf
        Mills_and_1000G_gold_standard.indels.hg38.vcf.gz

        [Protein] 
	
        This fold contains the reference cDNA and protein sequence of human:
        Homo_sapiens.GRCh38.pep.all.fa

6. Among the required software listed above, GATK, kallisto, picard, samtools, trimmomatic-0.36 were prepared in software directory, other software should be installed by user own due to complexity, please refer to the software links above.

7. Fill in the `config.yaml` file with your local path, make sure you have installed all above software and have downloaded reference data.You should be aware that the version of VEP library you use should match the references used (peptide and cDNA). E.g. in the example above used version/release 89 of GRCh38.


## Usage

You can use these two modes by:

        python PREP.py RE -i config.yaml


## Input Files

### Input Files
PREP accepts pair-end or single-end RNA sequencing as input. It could be in `.fastq.gz` or `.fastq` format. 
You should specify the right path to the sequencing file in `config.yaml` like:

    #your path to first RNA-seq fastq file
    tumor_rna_fastq_1: ~/ncbi/dbGaP-14145/sra/SRR2673065_1.fastq.gz
    #your path to second RNA-seq fastq file
    tumor_rna_fastq_2: ~/ncbi/dbGaP-14145/sra/SRR2673065_2.fastq.gz


## Setting parameters
User should set all the parameters in the configuration file `config.yaml` . The configuration file contains three parts of parameters:

* Input data parameters, including path of RNA sequencing data, output fold, run name.
* Software excutable path of opitype, vep, netMHCpan.

## Output Files 

The output files are the following: 
1.  final_neoantigen.tsv 

    The file is a TSV file with the extracted mutated peptides derived from RNA editing with a quantitative score measures the immunity of neoepitopes.


### Column explanation

The prediction output (final_neoantigen.tsv) for each peptide pair consists of the following columns:

| Column Name           | Description |
| -----------           | ----------- |
| Position              | Mutation position in genome. |
| HLA_type              | HLA allele name. |
| Gene                  | HUGO symbol name of mutatied gene. |
| Transcript_name       | Ensembl transcript ID |
| Mutation              | Necleotide change of mutated gene |
| AA_change             | Amino acid change annotated in VEP file. |
| WT_pep                | The extracted normal peptide. |
| WT_Binding_EL         | %Rank of prediction score for nomal peptides use NetMHCpan4.0 (defalut model). |
| WT_Binding_Rank       | %Rank of prediction score for nomal peptides use NetMHCpan4.0 (-ba model). |
| MT_pep                | The extracted mutant peptide. |
| MT_Binding_EL         | %Rank of prediction score for mutated peptides use NetMHCpan4.0(defalut model). |
| MT_Binding_Rank       | %Rank of prediction score for mutant peptides use NetMHCpan4.0 (-ba model). |
| DriverGene_Lable      | TRUE if the HUGO symbol is in the cosmic reference list, FALSE if it is not. |
| MT_Binding_level_des  | Binding level description of mutated peptide. |
| WT_Binding_level_des  | Binding level description of normal peptide. |
| Editing_ratio  | RNA ediitng level of the mutation. |
| Neo_score	| Immunogenicty score for RNA editing neoepitopes. |


## Contact   
 
1410782Chiz@tongji.edu.cn or qiliu@tongji.edu.cn

Biological and Medical Big data Mining Lab  
Tongji University    

## Algorithmic Flow Chart

![](PREP_flow_chart.jpg)
