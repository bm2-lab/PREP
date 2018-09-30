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
    

## Usage
1. Make sure `NetMHCpan` and `vep` are in your executable path.
2. Change `vep_cache` to your directory in `PREP.py` and Change `netchop_dir` to your directory in `bin/netCTLpan.py`.

For detailed usage, use:

    python path/to/PREP.py -h

## Test example

To run the provided test files with PREP the following command can be run: 

    python path/to/PREP.py -i test.vcf -o test -s test -a HLA-A02:01 -l 9 -e test.FPKM.txt -t test_editing_level.txt
    
## Contact   

Chi Zhou  
1410782Chiz@tongji.edu.cn 

Qi Liu  
qiliu@tongji.edu.cn

Biological and Medical Big data Mining Lab  
Tongji University    
Shanghai, China. 
