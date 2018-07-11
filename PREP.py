import pandas as pd
import numpy as np
import sys,getopt,os
import re
import os,sys,time
import multiprocessing
import shutil 
import math
from pyper import *
from sklearn import preprocessing
from sklearn.decomposition import PCA
from matplotlib import pyplot as plt
from sklearn.semi_supervised import label_propagation
import itertools
from scipy import linalg
import matplotlib as mpl
from sklearn import mixture
from sklearn.preprocessing import StandardScaler
from sklearn.neural_network import MLPClassifier
from imblearn.under_sampling import RandomUnderSampler
from imblearn.metrics import classification_report_imbalanced
from Bio.Blast import NCBIXML
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
from math import log, exp
from scipy import interp
from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import StratifiedKFold
import xgboost as xgb
from xgboost.sklearn import XGBClassifier
from sklearn import cross_validation, metrics
from sklearn.grid_search import GridSearchCV
import matplotlib.pylab as plt
from sklearn.model_selection import train_test_split
from Bio.Blast import NCBIXML
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
from math import log, exp
from sklearn.ensemble import RandomForestClassifier  
from sklearn.preprocessing import StandardScaler
from imblearn.over_sampling import SMOTE
from collections import Counter
from sklearn.cross_validation import cross_val_score
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
cmd_vep="vep -i " + input_vcf_file + " --cache --dir " + '/home/zhouchi/database/Annotation/vep_data_89_37/' + " --dir_cache " + '/home/zhouchi/database/Annotation/vep_data_89_37/' + " --force_overwrite --canonical --symbol -o STDOUT --offline | filter_vep --ontology --filter \"CANONICAL is YES and Consequence is missense_variant\" -o " + out_dir + '/' + sample_id +'_vep.txt' + " --force_overwrite"
print cmd_vep
os.system(cmd_vep)
cmd_edit="python bin/edit2fasta.py -i " + out_dir + "/" + sample_id +'_vep.txt' + " -o " + out_dir + " -s " + sample_id
print cmd_edit
os.system(cmd_edit)
hla_str=hla_allele	
cmd_netMHCpan="netMHCpan -a " + hla_str +  " -f " + out_dir + "/" + sample_id +'_edit.fasta' + " -l 9,10,11 -BA > " + out_dir + "/"  + sample_id +'_netmhc.txt'
print cmd_netMHCpan
os.system(cmd_netMHCpan)
cmd_parse="python bin/sm_netMHC_result_parse.py -i " + out_dir + "/"  + sample_id +'_netmhc.txt' + " -g " + out_dir + "/"  + sample_id +'_edit.fasta' + ' -t ' + out_dir + "/" + sample_id +'_vep.txt' + " -o " + out_dir + " -s " + sample_id + " -e " + exp_file_path + " -a 1 -b 2 -f 0 -l " + hla_str
print cmd_parse
os.system(cmd_parse)
cmd_netCTL="python bin/netCTLPAN.py -i " + out_dir + "/" + sample_id +"_final_neo_candidate.txt" + " -o " + out_dir + " -s " + sample_id
print cmd_netCTL
os.system(cmd_netCTL)

data_level=pd.read_table(editing_level_file,header = 0, sep='\t')
data_netclt=pd.read_table(out_dir + "/" + sample_id + "_netctl_concact.txt",header = 0, sep='\t')
data_out=pd.merge(data_netclt,data_level,left_on="#Position",right_on="POS",how="left")
del data_out["POS"]
#data_out.to_csv(out_dir + "/" + sample_id + "_neo.txt",header=1,index=0,sep='\t')

##some functions
###hydrophobicity
#preprocess
a=26
k=4.86936
M=1. #default concentration of mutant peptides
W=1. #default concentration of wildtype peptides

WEPS=0.0003
HYDROPHOBIC_RESIDUES="AILMFWYV"
WEIRD_RESIDUES="CGP"
hydro_score={"A":1.8,"C":2.5,"D":-3.5,"E":-3.5,"F":2.8,"G":-0.4,"H":-3.2,"I":4.5,"K":-3.9,"L":3.8,"M":1.9,"N":-3.5,"P":-1.6,"Q":-3.5,"R":-4.5,"S":-0.8,"T":-0.7,"V":4.2,"W":-0.9,"Y":-1.3,"X":0.0}



def get_iedb_seq(iedb_file):
	iedb_seq=[]
	for line in open(iedb_file):
		if line.startswith(">"):
			continue
		else:
			iedb_seq.append(line.strip())
	return iedb_seq
def GetNmerPositivePep(n,mhc_pos_file):
	pep_list=[]
	for line in open(mhc_pos_file):
		if line.startswith(">"):
			continue
		else:
			record=line.strip()
			if len(record)==n:
				pep_list.append(record)
	return pep_list


def GetNmerNegativePep(n,mhc_neg_file):
	pep_list=[]
	for line in open(mhc_neg_file):
		if line.startswith(">"):
			continue
		else:
			record=line.strip()
			if len(record)==n:
				pep_list.append(record)
	return pep_list

def hydro_vector(pep):
	hydrophobicity_vector=[]
	pep_list=list(pep)
	pep_len=len(pep_list)
	for pep in pep_list:
		hydrophobicity_vector.append(hydro_score[pep.upper()])

	return hydrophobicity_vector

def getXY(n,mhc_pos_file,mhc_neg_file):
	pos_pep=GetNmerPositivePep(n,mhc_pos_file)
	neg_pep=GetNmerNegativePep(n,mhc_neg_file)
	####get hydrophobicity list
	pos_pep_hydro_list=[hydro_vector(p) for p in pos_pep]#[0:len(neg_pep_all)]
	neg_pep_hydro_list=[hydro_vector(p) for p in neg_pep]
	#print len(pos_pep_hydro_list)
	#print len(neg_pep_hydro_list)
	######transform into array
	#pos_pep_hydro_array=np.array(pos_pep_hydro_list)
	#neg_pep_hydro_array=np.array(neg_pep_hydro_list)
	X=pos_pep_hydro_list+neg_pep_hydro_list
	y=[1]*len(pos_pep_hydro_list)+[0]*len(neg_pep_hydro_list)
	X_array=np.asarray(X)
	y_array=np.asarray(y)
	return X_array,y_array

###get model
def get_hydro_model(mhc_pos_file,mhc_neg_file):
	xgb_9 = XGBClassifier(
	 learning_rate =0.1,
	 n_estimators=208,
	 max_depth=5,
	 min_child_weight=1,
	 gamma=0,
	 subsample=0.8,
	 colsample_bytree=0.8,
	 objective= 'binary:logistic',
	 n_jobs=4,
	 scale_pos_weight=1,
	 random_state=27)
	#[207]  train-auc:0.95222+0.000337282   test-auc:0.852987+0.00777271
	#AUC Score (Train): 0.941597
	#AUC Score (Test): 0.768922

	xgb_10 = XGBClassifier(
	 learning_rate =0.1,
	 n_estimators=165,
	 max_depth=16,
	 min_child_weight=1,
	 gamma=0,
	 subsample=0.8,
	 colsample_bytree=0.8,
	 objective= 'binary:logistic',
	 n_jobs=4,
	 scale_pos_weight=1,
	 random_state=27)
	#[164]  train-auc:1+0   test-auc:0.900267+0.00955733
	#AUC Score (Train): 1.000000
	#AUC Score (Test): 0.802474
	xgb_11 = XGBClassifier(
	 learning_rate =0.1,
	 n_estimators=123,
	 max_depth=5,
	 min_child_weight=1,
	 gamma=0.1,
	 subsample=0.8,
	 colsample_bytree=0.8,
	 objective= 'binary:logistic',
	 n_jobs=4,
	 scale_pos_weight=1,
	 random_state=27)
	X_array_9,y_array_9=getXY(9,mhc_pos_file,mhc_neg_file)
	hy_xgb_9=xgb_9.fit(X_array_9,y_array_9)
	X_array_10,y_array_10=getXY(10,mhc_pos_file,mhc_neg_file)
	hy_xgb_10=xgb_10.fit(X_array_10,y_array_10)
	X_array_11,y_array_11=getXY(11,mhc_pos_file,mhc_neg_file)
	hy_xgb_11=xgb_11.fit(X_array_11,y_array_11)
	return hy_xgb_9,hy_xgb_10,hy_xgb_11

######T cell recognition
def logSum(v):
	ma=max(v)
	return log(sum(map(lambda x: exp(x-ma),v)))+ma



def calculate_R(neo_seq,iedb_seq):
	align_score=[]
	#i=0
	for seq in iedb_seq:
		aln_score=aligner(neo_seq,seq)
		#i=i+1
		#print i
		#print aln_score
		if aln_score!=[]:
			localds_core=max([line[2] for line in aln_score])
			align_score.append(localds_core)
	#print align_score
	#print k,a
	bindingEnergies=map(lambda x: -k*(a-x),align_score)
	#print bindingEnergies
	lZk=logSum(bindingEnergies+[0])
	lGb=logSum(bindingEnergies)
	R=exp(lGb-lZk)
	return R

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
iedb_seq=get_iedb_seq('/home/zhouchi/database/Annotation/protein/iedb.fasta')
hy_xgb_9,hy_xgb_10,hy_xgb_11=get_hydro_model('/home/zhouchi/database/Annotation/protein/all_postive.txt','/home/zhouchi/database/Annotation/protein/all_negative.txt')

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
	#r=calculate_R(MT_peptide[i],iedb_seq)
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
data_out_score.to_csv(out_dir + "/" + sample_id + "_neo_rank.txt",header=1,index=0,sep='\t')
