# PREP: Prioritizing of RNA Editing-based Peptides

RNA editing is a source of transcriptomic diversity, mainly in non-coding regions, and is altered in cancer. Recent studies demonstrated that A-to-I RNA editing events are manifested at the proteomic level and contribute to protein heterogeneity in cancer. Given somatic RNA-editing mutation as input, PREP identify and evaluate the potential immunogenicity of RNA editing based peptides. Detailed information please refer to citation.

#### Authors:
Chi Zhou and Qi Liu

#### Citation:
A-to-I RNA Editing Generates Potential Neoantigens in Cancer, Submitted, 2018.

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


### Install from source
1. Install all software and python package listed above.

2. Download or clone the PREP repository to your local system:

        git clone https://github.com/bm2-lab/PREP.git

2. Change `NetMHCpan`, `vep` and `vep_cache` to your directory in `config.yaml`. Please refer to [user manual](/doc/PREP_User_Manual.md) for a detailed description.



## Usage

To run PREP, first fill file path of you rna-seq FASTQ file in `confg.yaml`, then use:

    python PREP.py RE -i config.yaml
  
## User Manual 
For detailed information about usage, input and output files, test examples and data
preparation please refer to the [PREP User Manual](/doc/PREP_User_Manual.md)

## Contact   

Chi Zhou  
1410782Chiz@tongji.edu.cn 

Qi Liu  
qiliu@tongji.edu.cn

Biological and Medical Big data Mining Lab  
Tongji University    
Shanghai, China. 
