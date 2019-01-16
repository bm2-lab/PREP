import sys, os

input_vcf=sys.argv[1]
filter_vcf=sys.argv[2]


normal_list=[]
with open('bin/normal_editing_set.txt','r') as f_i:
	for line in f_i:
		normal_list.append(line.strip())

with open(input_vcf,'r') as fin, open(filter_vcf,'w') as fout:
	for line in fin:
		rec=line.strip().split('\t')
		c_c=rec[0]+':'+rec[1]
		if c_c not in normal_list:
			fout.write(line)