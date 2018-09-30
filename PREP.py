import pandas as pd
import numpy as np
import sys,getopt,os
import re
import math
from sklearn import preprocessing
from sklearn.preprocessing import StandardScaler
from sklearn.neural_network import MLPClassifier
from Bio.Blast import NCBIXML
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
from math import log, exp
import xgboost as xgb
from xgboost.sklearn import XGBClassifier
from sklearn.ensemble import RandomForestClassifier  
from sklearn.externals import joblib

#####prepare fasta format input file for netMHC######
opts,args=getopt.getopt(sys.argv[1:],"hi:o:s:a:l:e:t:",["input_vcf_file","out_dir","sample_id","hla_allele","pep_len","exp_file_path","editing_level_file"])
input_vcf_file =""
out_dir=""
sample_id=""
hla_allele=""
pep_len="9"
exp_file_path="no_exp"
editing_level_file=""
USAGE='''
This script convert RNA editing site to candidate neoantigen
usage: python RNAEditing.py -i <input_vcf_file> -o <outdir> -s <sample_id> -a <hla_allele> -l <pep_len> -e <exp_file_path> -t <editing_level_file>
required argument:
		-i | --input_vcf_file : input file containing RNA editing sites
		-o | --out_dir : output directory
		-s | --sample_id : sample id
		-a | --hla_allele : hla allele, separate by comma, eg: HLA-A02:01,HLA-A03:02
		-l | --pep_len : candidate neo-peptide length, separate by comma,eg: 9,10,11. default: 9
		-e | --exp_file_path : expression file
		-t | --editing_level_file : RNA editing level file
'''
for opt,value in opts:
	if opt =="h":
		print USAGE
		sys.exit(2)
	elif opt in ("-i","--input_vcf_file"):
		input_vcf_file=value
	elif opt in ("-o","--out_dir"):
		out_dir =value
	elif opt in ("-s","--sample_id"):
		sample_id =value  
	elif opt in ("-a","--hla_allele"):
		hla_allele =value 
	elif opt in ("-l","--pep_len"):
		pep_len =value 
	elif opt in ("-e","--exp_file_path"):
		exp_file_path =value
	elif opt in ("-t","--editing_level_file"):
		editing_level_file =value	
#print coverage
if (input_vcf_file =="" or out_dir =="" or sample_id=="" or hla_allele==""):
	print USAGE
	sys.exit(2)

##########change software path to your own directory###########
#vep_cache='your/path/to/vep cache data'#'/home/zhouchi/database/Annotation/vep_data_89_37/'
#vep_path='your/path/to/vep'
#netMHCpan_path='your/path/to/netMHCpan'

vep_cache='/home/zhouchi/database/Annotation/vep_data_89_37/'
vep_path='vep'
netMHCpan_path='netMHCpan'
netchop_path='/home/zhouchi/software/netchop'


cmd_vep=vep_path + " -i " + input_vcf_file + " --cache --dir " + vep_cache + " --dir_cache " + vep_cache + " --force_overwrite --canonical --symbol -o STDOUT --offline | filter_vep --ontology --filter \"CANONICAL is YES and Consequence is missense_variant\" -o " + out_dir + '/' + sample_id +'_vep.tsv' + " --force_overwrite"
print "Mutation annotation using VEP..."
os.system(cmd_vep)
cmd_edit="python bin/edit2fasta.py -i " + out_dir + "/" + sample_id +'_vep.tsv' + " -o " + out_dir + " -s " + sample_id
print "Convert mutation to fasta..."
os.system(cmd_edit)
hla_str=hla_allele	
cmd_netMHCpan=netMHCpan_path + " -a " + hla_str +  " -f " + out_dir + "/" + sample_id +'_edit.fasta' + " -l " + pep_len + " -BA > " + out_dir + "/"  + sample_id +'_netmhc.tsv'
print "Neoantigen identification using netMHCpan..."
os.system(cmd_netMHCpan)
cmd_parse="python bin/sm_netMHC_result_parse.py -i " + out_dir + "/"  + sample_id +'_netmhc.tsv' + " -g " + out_dir + "/"  + sample_id +'_edit.fasta' + ' -t ' + out_dir + "/" + sample_id +'_vep.tsv' + " -o " + out_dir + " -s " + sample_id + " -e " + exp_file_path + " -a 1 -b 2 -f 0 -l " + hla_str
print "Neoantigen extraction..."
os.system(cmd_parse)
cmd_netCTL="python bin/netCTLPAN.py -i " + out_dir + "/" + sample_id +"_final_neo_candidate.tsv" + " -d " + "data/DriveGene.tsv " + " -o " + out_dir + " -s " + sample_id + " -n " + netchop_path
print "Neoantigen annotation and ranking by immunogenicity score..."
os.system(cmd_netCTL)

data_level=pd.read_table(editing_level_file,header = 0, sep='\t')
data_netclt=pd.read_table(out_dir + "/" + sample_id + "_netctl_concact.tsv",header = 0, sep='\t')
data_out=pd.merge(data_netclt,data_level,left_on="#Position",right_on="POS",how="left")
del data_out["POS"]
#data_out.to_csv(out_dir + "/" + sample_id + "_neo.txt",header=1,index=0,sep='\t')

##some functions
###hydrophobicity
#preprocess

HYDROPHOBIC_RESIDUES="AILMFWYV"
WEIRD_RESIDUES="CGP"
hydro_score={"A":1.8,"C":2.5,"D":-3.5,"E":-3.5,"F":2.8,"G":-0.4,"H":-3.2,"I":4.5,"K":-3.9,"L":3.8,"M":1.9,"N":-3.5,"P":-1.6,"Q":-3.5,"R":-4.5,"S":-0.8,"T":-0.7,"V":4.2,"W":-0.9,"Y":-1.3,"X":0.0}


def hydro_vector(pep):
	hydrophobicity_vector=[]
	pep_list=list(pep)
	pep_len=len(pep_list)
	for pep in pep_list:
		hydrophobicity_vector.append(hydro_score[pep.upper()])

	return hydrophobicity_vector



############similarity############
def aligner(seq1,seq2):
	matrix = matlist.blosum62
	gap_open = -11
	gap_extend = -1
	aln = pairwise2.align.localds(seq1.upper(), seq2.upper(), matrix,gap_open,gap_extend)
	#print aln
	return aln

def cal_similarity_per(mut_seq,normal_seq):
	score_pair=aligner(mut_seq,normal_seq)[0][2]
	score_self=aligner(mut_seq,mut_seq)[0][2]
	per_similarity=score_pair/score_self
	return per_similarity

hy_xgb_9=joblib.load("model/cf_hy_9_model.m")
hy_xgb_10=joblib.load("model/cf_hy_10_model.m")
hy_xgb_11=joblib.load("model/cf_hy_11_model.m")
fpkm=data_out.FPKM
MT_peptide=data_out.MT_pep
WT_peptide=data_out.WT_pep
MT_rank=data_out.MT_Binding_Aff
WT_rank=data_out.WT_Binding_Aff
MT_editing_level=data_out.editing_level
f_fpkm=lambda x:math.tanh(x)
fpkm_score=list(data_out.FPKM.apply(f_fpkm))
f_rank_wt=lambda x:1-(1/(1+math.pow(math.e,5*(float(x)-2))))/2
f_rank_mt=lambda x:1/(1+math.pow(math.e,5*(float(x)-2)))
mt_rank_score=list(data_out.MT_Binding_Aff.apply(f_rank_mt))
wt_rank_score=list(data_out.WT_Binding_Aff.apply(f_rank_wt))
immuno_score=[]
for i in range(len(MT_peptide)):
	s=1.0-cal_similarity_per(MT_peptide[i],WT_peptide[i])
	e=data_out.editing_level[i]
	pep_len=len(MT_peptide[i])
	if pep_len==9:
		h=hy_xgb_9.predict_proba(np.array(hydro_vector(MT_peptide[i])).reshape((1,9)))[:,1][0]
	elif pep_len==10:
		h=hy_xgb_10.predict_proba(np.array(hydro_vector(MT_peptide[i])).reshape((1,10)))[:,1][0]
	else:
		h=hy_xgb_11.predict_proba(np.array(hydro_vector(MT_peptide[i])).reshape((1,11)))[:,1][0]
	score_i=e*h*s*mt_rank_score[i]*wt_rank_score[i]*fpkm_score[i]
	immuno_score.append(score_i)
data_out["immuno_score"]=immuno_score
data_out_score=data_out.sort_values(["immuno_score"],ascending=[0])
data_out_score.to_csv(out_dir + "/" + sample_id + "_neo_rank.tsv",header=1,index=0,sep='\t')
print "All finished, check neo_rank.tsv for results!"