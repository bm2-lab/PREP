# PREP: Prioritizing of RNA Editing-based Peptides

Given somatic RNA-editing mutation as input, PREP identify and evaluate the potential immunogenicity of RNA editing based peptides. Detailed information please refer to citation.

#### Citation:
Somatic A-to-I RNA Editing Generates Potential Neoantigens in Cancer, Submitted, 2018.

## Dependencies  

#### Required software:
* [NetMHCpan 4.0](http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?netMHCpan)
* [Variant Effect Predictor](https://github.com/Ensembl/ensembl-vep)
* [NetChop](http://www.cbs.dtu.dk/services/NetChop/)

#### Python package
    xgboost
    biopython
    scikit-learn
    pandas
    numpy
    
### Install via Docker
Docker image of PREP and a detailed description usage is at https://hub.docker.com/r/bm2lab/prep/.


### Install from source
1. Install all software and python package listed above.

2. Download or clone the PREP repository to your local system:

        git clone https://github.com/bm2-lab/PREP.git

2. Change `NetMHCpan`, `vep`, `netchop` and `vep_cache` to your directory in `PREP.py`.



## Usage

For detailed usage, use:

    python path/to/PREP.py -h

## Test example

To run the provided test files with PREP the following command can be run: 

    python PREP.py -i example/example.vcf -o example -s example -a HLA-A02:01 -l 9 -e example_FPKM.txt -t example_editing_level.txt
    
## Contact   

Chi Zhou  
1410782Chiz@tongji.edu.cn 

Qi Liu  
qiliu@tongji.edu.cn

Biological and Medical Big data Mining Lab  
Tongji University    
Shanghai, China. 

## Neoantigen presentation illustration and algorithmic flow chart

![](/doc/PREP_flow_chart.jpg)
