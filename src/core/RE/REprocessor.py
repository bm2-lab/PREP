import os,sys,time
import multiprocessing
import shutil 
import subprocess
import pandas as pd
import math
from pyper import *
import numpy as np
import math
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
import pandas as pd
import math	
import numpy as np
from scipy import interp
from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import StratifiedKFold
import xgboost as xgb
from xgboost.sklearn import XGBClassifier
from sklearn import metrics
from sklearn.model_selection import GridSearchCV
import matplotlib.pylab as plt
from sklearn.model_selection import train_test_split
from Bio.Blast import NCBIXML
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
from math import log, exp
import subprocess
from sklearn.ensemble import RandomForestClassifier  
from sklearn.preprocessing import StandardScaler
from imblearn.over_sampling import SMOTE
from collections import Counter
from sklearn.model_selection import cross_val_score
from sklearn.externals import joblib

a=26
k=4.86936
M=1. #default concentration of mutant peptides
W=1. #default concentration of wildtype peptides

WEPS=0.0003
HYDROPHOBIC_RESIDUES="AILMFWYV"
WEIRD_RESIDUES="CGP"
def get_iedb_seq(iedb_file):
	iedb_seq=[]
	for line in open(iedb_file):
		if line.startswith(">"):
			continue
		else:
			iedb_seq.append(line.strip())
	return iedb_seq
def read_trimmomatic_PE(raw_fastq_path_first,raw_fastq_path_second,trimmomatic_path,adapter_path,fastq_prefix,logfile_fold,CPU):
	cmd_trimmomatic="java -jar " + trimmomatic_path + " PE -phred33 -threads " + str(CPU) + ' ' + raw_fastq_path_first + ' ' + raw_fastq_path_second  + ' -baseout ' + fastq_prefix + " ILLUMINACLIP:" + adapter_path + ':2:30:10' + ' LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 > ' + logfile_fold + '/' + 'trimmomatic.log' + ' 2>&1'
	#print cmd_trimmomatic
	os.system(cmd_trimmomatic)

def read_trimmomatic_SE(raw_fastq_path_first,trimmomatic_path,adapter_path,fastq_prefix,logfile_fold,CPU):
	cmd_trimmomatic="java -jar " + trimmomatic_path + " SE -phred33 -threads " + str(CPU) + ' ' + raw_fastq_path_first + ' -baseout ' + fastq_prefix + " ILLUMINACLIP:" + adapter_path + ':2:30:10' + ' LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 > ' + logfile_fold + '/' + 'trimmomatic.log' + ' 2>&1'
	#print cmd_trimmomatic
	os.system(cmd_trimmomatic)


def hlatyping_pe(raw_fastq_path_first,raw_fastq_path_second,opitype_fold,opitype_out_fold,opitype_ext,prefix):
	cmd_hla = 'python ' + opitype_fold + ' -i ' + raw_fastq_path_first + ' ' + raw_fastq_path_second + ' --rna -o ' + opitype_out_fold
	#print cmd_hla
	os.system(cmd_hla)
	result_dir=os.listdir(opitype_out_fold)
	#print result_dir[0]
	hla_result_path=opitype_out_fold+'/'+result_dir[0]+'/'+result_dir[0]+'_result.tsv'
	#print hla_result_path
	cmd_hla_ext = 'python ' + opitype_ext + ' -i ' + hla_result_path + ' -o ' + opitype_out_fold + ' -s ' + prefix
	#print cmd_hla_ext
	os.system(cmd_hla_ext)
	print 'hla type process done.'

def hlatyping_se(raw_fastq_path_first,opitype_fold,opitype_out_fold,opitype_ext,prefix):
	cmd_hla = 'python ' + opitype_fold + ' -i ' + raw_fastq_path_first + ' --rna -o ' + opitype_out_fold
	#print cmd_hla
	os.system(cmd_hla)
	result_dir=os.listdir(opitype_out_fold)
	#print result_dir[0]
	hla_result_path=opitype_out_fold+'/'+result_dir[0]+'/'+result_dir[0]+'_result.tsv'
	#print hla_result_path
	cmd_hla_ext = 'python ' + opitype_ext + ' -i ' + hla_result_path + ' -o ' + opitype_out_fold + ' -s ' + prefix
	#print cmd_hla_ext
	os.system(cmd_hla_ext)
	print 'hla type process done.'





def mapping_PE(fastq_1_path,fastq_2_path,CPU,alignment_out_fold,prefix,star_path,star_index,stringtie_path,gtf_path,picard_path,gatk_path,reference,indels,rnaeditor_path):
	cmd_star="{} --twopassMode Basic --genomeDir {} --runThreadN {} --outSAMtype BAM SortedByCoordinate --twopass1readsN -1 --sjdbOverhang 99 --readFilesIn {} {} --outSAMattrRGline ID:RG_{} SM:{} PL:ILLUMINA --limitBAMsortRAM 60000000000 --outFileNamePrefix {}".format(star_path,star_index,CPU,fastq_1_path,fastq_2_path,prefix,prefix,alignment_out_fold+'/'+prefix)
	print cmd_star
	#os.system(cmd_star)
	cmd_stringtie="{} -G {} -A {}_exp.gtf -o {}.gtf {}Aligned.sortedByCoord.out.bam".format(stringtie_path,gtf_path,alignment_out_fold+'/'+prefix,alignment_out_fold+'/'+prefix,alignment_out_fold+'/'+prefix)
	print cmd_stringtie
	#os.system(cmd_stringtie)
	cmd_mkdup="java -jar {} MarkDuplicates I={}Aligned.sortedByCoord.out.bam O={}_markdup.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M={}.metrics".format(picard_path,alignment_out_fold+'/'+prefix,alignment_out_fold+'/'+prefix,alignment_out_fold+'/'+prefix)
	print cmd_mkdup
	#os.system(cmd_mkdup)
	cmd_splitjunc="java -jar {} -T SplitNCigarReads -R {} -I {}_markdup.bam -o {}_splitted.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS".format(gatk_path,reference,alignment_out_fold+'/'+prefix,alignment_out_fold+'/'+prefix)
	print cmd_splitjunc
	#os.system(cmd_splitjunc)
	cmd_RTC="java -jar {} -T RealignerTargetCreator -R {} -I {}_splitted.bam --known {} -o {}.intervalListFromRTC.intervals".format(gatk_path,reference,alignment_out_fold+'/'+prefix,indels,alignment_out_fold+'/'+prefix)
	print cmd_RTC
	#os.system(cmd_RTC)		
	cmd_indelrealign="java -jar {} -T IndelRealigner -R {} -I {}_splitted.bam -o {}_realigned.bam -known {} -targetIntervals {}.intervalListFromRTC.intervals".format(gatk_path,reference,alignment_out_fold+'/'+prefix,alignment_out_fold+'/'+prefix,indels,alignment_out_fold+'/'+prefix)
	print cmd_indelrealign
	#os.system(cmd_indelrealign)
	cmd_rename_bam="mv {}_realigned.bam {}.noDup.realigned.recalibrated.bam".format(alignment_out_fold+'/'+prefix,alignment_out_fold+'/'+prefix)
	cmd_rename_bai="mv {}_realigned.bai {}.noDup.realigned.recalibrated.bai".format(alignment_out_fold+'/'+prefix,alignment_out_fold+'/'+prefix)
	print cmd_rename_bam
	#os.system(cmd_rename_bam)
	#os.system(cmd_rename_bai)
	cmd_chcon="sed '11c output = {}' {}configuration.txt > {}configuration_new.txt".format(alignment_out_fold+'/'+prefix,rnaeditor_path,rnaeditor_path)
	print cmd_chcon
	os.system(cmd_chcon)
	cmd_rnaeditor="python {}RNAEditor.py -i {}.noDup.realigned.recalibrated.bam -c {}configuration_new.txt".format(rnaeditor_path,alignment_out_fold+'/'+prefix,rnaeditor_path)
	print cmd_rnaeditor
	os.system(cmd_rnaeditor)

def mapping_SE(fastq_1_path,CPU,alignment_out_fold,prefix,star_path,star_index,stringtie_path,gtf_path,picard_path,gatk_path,reference,indels,rnaeditor_path):
	cmd_star="{} --twopassMode Basic --genomeDir {} --runThreadN {} --outSAMtype BAM SortedByCoordinate --twopass1readsN -1 --sjdbOverhang 99 --readFilesIn {} {} --outSAMattrRGline ID:RG_{} SM:{} PL:ILLUMINA --limitBAMsortRAM 60000000000 --outFileNamePrefix {}".format(star_path,star_index,CPU,fastq_1_path,prefix,prefix,alignment_out_fold+'/'+prefix)
	print cmd_star
	os.system(cmd_star)
	cmd_stringtie="{} -G {} -A {}_exp.gtf -o {}.gtf {}Aligned.sortedByCoord.out.bam".format(stringtie_path,gtf_path,alignment_out_fold+'/'+prefix,alignment_out_fold+'/'+prefix,alignment_out_fold+'/'+prefix)
	print cmd_stringtie
	os.system(cmd_stringtie)
	cmd_mkdup="java -jar {} MarkDuplicates I={}Aligned.sortedByCoord.out.bam O={}_markdup.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M={}.metrics".format(picard_path,alignment_out_fold+'/'+prefix,alignment_out_fold+'/'+prefix,alignment_out_fold+'/'+prefix)
	print cmd_mkdup
	#os.system(cmd_mkdup)
	cmd_splitjunc="java -jar {} -T SplitNCigarReads -R {} -I {}_markdup.bam -o {}_splitted.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS".format(gatk_path,reference,alignment_out_fold+'/'+prefix,alignment_out_fold+'/'+prefix)
	print cmd_splitjunc
	#os.system(cmd_splitjunc)
	cmd_RTC="java -jar {} -T RealignerTargetCreator -R {} -I {}_splitted.bam --known {} -o {}.intervalListFromRTC.intervals".format(gatk_path,reference,alignment_out_fold+'/'+prefix,indels,alignment_out_fold+'/'+prefix)
	print cmd_RTC
	#os.system(cmd_RTC)		
	cmd_indelrealign="java -jar {} -T IndelRealigner -R {} -I {}_splitted.bam -o {}_realigned.bam -known {} -targetIntervals {}.intervalListFromRTC.intervals".format(gatk_path,reference,alignment_out_fold+'/'+prefix,alignment_out_fold+'/'+prefix,indels,alignment_out_fold+'/'+prefix)
	print cmd_indelrealign
	#os.system(cmd_indelrealign)
	cmd_rename_bam="mv {}_realigned.bam {}.noDup.realigned.recalibrated.bam".format(alignment_out_fold+'/'+prefix,alignment_out_fold+'/'+prefix)
	cmd_rename_bai="mv {}_realigned.bai {}.noDup.realigned.recalibrated.bai".format(alignment_out_fold+'/'+prefix,alignment_out_fold+'/'+prefix)
	print cmd_rename_bam
	#os.system(cmd_rename_bam)
	#os.system(cmd_rename_bai)
	#cmd_chcon="sed '11c output = {}' {}configuration.txt > {}configuration_new.txt".format(alignment_out_fold+'/'+prefix,rnaeditor_path,rnaeditor_path)
	#print cmd_chcon
	#os.system(cmd_chcon)
	cmd_rnaeditor="python {}RNAEditor.py -i {}.noDup.realigned.recalibrated.bam -c {}configuration_new.txt".format(rnaeditor_path,alignment_out_fold+'/'+prefix,rnaeditor_path)
	print cmd_rnaeditor
	#os.system(cmd_rnaeditor)


def neo_detection(alignment_out_fold,prefix,vep_path,vep_cache,human_peptide_path,netmhcpan_path,hla_str,netchop_path):
	cmd_coding="grep -w \"coding-exon\" {}.editingSites.vcf | awk -F \'\\t\' '{{print \"chr\"$0}}' - > {}.CodingeditingSites_chr.vcf".format(alignment_out_fold+'/'+prefix,alignment_out_fold+'/'+prefix)
	print cmd_coding
	os.system(cmd_coding)
	cmd_filter_normal="python bin/editing_site_filtering.py {}.CodingeditingSites_chr.vcf {}.CodingeditingSites_chr_filter.vcf".format(alignment_out_fold+'/'+prefix,alignment_out_fold+'/'+prefix)
	print cmd_filter_normal
	os.system(cmd_filter_normal)
	cmd_vep="{}vep -i {}.CodingeditingSites_chr_filter.vcf --cache --dir {} --dir_cache {} --force_overwrite --canonical --symbol -o STDOUT --offline | {}filter_vep --ontology --filter 'Consequence is missense_variant' -o {}.vep.txt --force_overwrite".format(vep_path,alignment_out_fold+'/'+prefix,vep_cache,vep_cache,vep_path,alignment_out_fold+'/'+prefix)
	print cmd_vep
	os.system(cmd_vep)
	cmd_gene_fasta="python bin/snv2fasta.py -i {}.vep.txt -o {} -s {} -p {}".format(alignment_out_fold+'/'+prefix,alignment_out_fold,prefix,human_peptide_path)
	print cmd_gene_fasta
	os.system(cmd_gene_fasta)
	cmd_netmhcpan="{} -f {}_snv.fasta -a {} -l 9,10,11 -BA > {}_nemhc.txt".format(netmhcpan_path,alignment_out_fold+'/'+prefix,hla_str,alignment_out_fold+'/'+prefix)
	print cmd_netmhcpan
	os.system(cmd_netmhcpan)
	cmd_netmhcpan_parse="python bin/sm_netMHC_result_parse.py -i {}_nemhc.txt -g {}_snv.fasta -o {} -s {} -l {} -e {}_exp.gtf".format(alignment_out_fold+'/'+prefix,alignment_out_fold+'/'+prefix,alignment_out_fold,prefix,hla_str,alignment_out_fold+'/'+prefix)
	print cmd_netmhcpan_parse
	os.system(cmd_netmhcpan_parse)
	cmd_netchop="python bin/netCTLPAN.py -i {}_final_neo_candidate.tsv -o {} -d software/DriveGene.tsv -s {} -n {}".format(alignment_out_fold+'/'+prefix,alignment_out_fold,prefix,netchop_path)
	print cmd_netchop
	os.system(cmd_netchop)


def hydro_vector(pep):
	hydro_score={"A":1.8,"C":2.5,"D":-3.5,"E":-3.5,"F":2.8,"G":-0.4,"H":-3.2,"I":4.5,"K":-3.9,"L":3.8,"M":1.9,"N":-3.5,"P":-1.6,"Q":-3.5,"R":-4.5,"S":-0.8,"T":-0.7,"V":4.2,"W":-0.9,"Y":-1.3}
	hydrophobicity_vector=[]
	pep_list=list(pep)
	pep_len=len(pep_list)
	for pep in pep_list:
		hydrophobicity_vector.append(hydro_score[pep.upper()])
	return hydrophobicity_vector



def aligner(seq1,seq2):
	matrix = matlist.blosum62
	gap_open = -11
	gap_extend = -1
	aln = pairwise2.align.localds(seq1.upper(), seq2.upper(), matrix,gap_open,gap_extend)
	#print aln
	return aln



############similarity############

def cal_similarity_per(mut_seq,normal_seq):
	score_pair=aligner(mut_seq,normal_seq)[0][2]
	score_self=aligner(mut_seq,mut_seq)[0][2]
	per_similarity=score_pair/score_self
	return per_similarity


def neo_score(editing_gvf_file,neo_file,final_neo_file,cf_hy_model_9,cf_hy_model_10,cf_hy_model_11):
	hy_xgb_9=joblib.load(cf_hy_model_9)
	hy_xgb_10=joblib.load(cf_hy_model_10)
	hy_xgb_11=joblib.load(cf_hy_model_11)
	editing_ratio_dic={}
	with open(editing_gvf_file,'r') as fin:
		fin.readline()
		for line in fin.readlines():
			rec=line.strip().split('\t')
			cor="chr"+rec[3]+":"+rec[7]
			editng_ratio=rec[17]
			if rec[2].split(',')[0]=="coding-exon":
				editing_ratio_dic[cor]=editng_ratio	
	data_neo=pd.read_table(neo_file,sep='\t',header=0)
	editing_list=[]
	for line in data_neo["#Position"]:
		editing_list.append(editing_ratio_dic[line])
	data_neo["editing_ratio"]=editing_list
	tpm=data_neo.TPM
	MT_peptide=data_neo.MT_pep
	WT_peptide=data_neo.WT_pep
	MT_rank=data_neo.MT_Binding_Aff
	WT_rank=data_neo.WT_Binding_Aff
	MT_editing_level=data_neo.editing_ratio
	f_fpkm=lambda x:math.tanh(x)
	fpkm_score=list(data_neo.TPM.apply(f_fpkm))
	f_rank_wt=lambda x:1-(1/(1+math.pow(math.e,5*(float(x)-2))))/2
	f_rank_mt=lambda x:1/(1+math.pow(math.e,5*(float(x)-2)))
	mt_rank_score=list(data_neo.MT_Binding_Aff.apply(f_rank_mt))
	wt_rank_score=list(data_neo.WT_Binding_Aff.apply(f_rank_wt))
	hydrophobicity=[]
	dissimilarity=[]
	score_list=[]
	for i in range(len(MT_peptide)):
		s=1.0-cal_similarity_per(MT_peptide[i],WT_peptide[i])
		e=data_neo.editing_ratio[i]
		pep_len=len(MT_peptide[i])
		if pep_len==9:
			h=hy_xgb_9.predict_proba(np.array(hydro_vector(MT_peptide[i])).reshape((1,9)))[:,1][0]
		elif pep_len==10:
			h=hy_xgb_10.predict_proba(np.array(hydro_vector(MT_peptide[i])).reshape((1,10)))[:,1][0]
		else:
			h=hy_xgb_11.predict_proba(np.array(hydro_vector(MT_peptide[i])).reshape((1,11)))[:,1][0]
		dissimilarity.append(s)
		hydrophobicity.append(h)
		#print e,h,s,mt_rank_score[i],wt_rank_score[i],fpkm_score[i]
		score_i=float(e)*float(h)*float(s)*float(mt_rank_score[i])*float(wt_rank_score[i])*float(fpkm_score[i])
		#print score_i
		score_list.append(score_i)
	data_neo["neo_score"]=score_list
	data_neo.to_csv(final_neo_file,sep='\t',header=1,index=0)










