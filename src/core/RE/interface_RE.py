import os
from REprocessor import *
import shutil
import yaml
import time
def RNAEditing(opts):
	config_file=opts.Config_file
	f=open(config_file)
	config_list=yaml.load(f)
	#######read and parse parameter
	print "Read and parse parameters..."
	output_fold=config_list["output_fold"]
	itunes_bin_path="bin"
	opitype_fold=config_list["opitype_fold"]
	opitype_out_fold=output_fold + '/' + 'hlatyping'
	opitype_ext=itunes_bin_path+'/optitype_ext.py'
	prefix=config_list["sample_name"]
	CPU=config_list["thread_number"]
	vcftools_path="software/vcftools"
	REFERENCE="database/GRCH38/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
	star_index_path="database/GRCH38/STAR_index"
	alignment_out_fold=output_fold + '/' + 'alignments'
	star_path="software/STAR"
	samtools_path="software/samtools"
	java_picard_path="software/picard.jar"
	GATK_path="software/GATK/GenomeAnalysisTK.jar"
	somatic_mutation_fold=output_fold + '/' + 'somatic_mutation'
	vep_cache=config_list["vep_cache_path"]
	vep_path=config_list["vep_path"]
	clean_fastq_fold=output_fold + '/' + 'clean_fastq'
	netmhc_out_fold=output_fold + '/' + 'netmhc'
	logfile_out_fold=output_fold + '/' + 'logfile'
	human_peptide_path="database/GRCH38/Homo_sapiens.GRCh38.pep.all.fa"
	netctl_out_fold=output_fold + '/' + 'netctl'
	netMHCpan_path=config_list["netMHCpan_path"]
	snv_fasta_file=netmhc_out_fold+'/'+prefix+'_snv.fasta'
	snv_netmhc_out_file=netmhc_out_fold+'/'+prefix+'_snv_netmhc.tsv'
	rna_fastq_1_path=config_list["tumor_rna_fastq_1"]
	rna_fastq_2_path=config_list["tumor_rna_fastq_2"]
	trimmomatic_path="software/trimmomatic-0.36.jar"
	adapter_path_PE="software/TruSeq3-PE.fa"
	adapter_path_SE="software/TruSeq3-SE.fa"
	tumor_fastq_prefix=clean_fastq_fold + '/' + prefix + "_tumor_clean.fq"
	tumor_fastq_clean_1=clean_fastq_fold + '/' + prefix + "_tumor_clean_1P.fq"
	tumor_fastq_clean_2=clean_fastq_fold + '/' + prefix + "_tumor_clean_2P.fq"
	cf_hy_model_9="train_model/cf_hy_9_model.m"
	cf_hy_model_10="train_model/cf_hy_10_model.m"
	cf_hy_model_11="train_model/cf_hy_11_model.m"
	stringtie_path="software/stringtie"
	gtf_path="database/GRCH38/Homo_sapiens.GRCh38.83.gtf"
	indels="database/GRCH38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
	netchop_path="software/netchop"
	editing_gvf_file=alignment_out_fold+"/"+prefix+".editingSites.gvf"
	neo_file=alignment_out_fold+"/"+prefix+"_netctl_concact.tsv"
	final_neo_file=alignment_out_fold+"/"+prefix+"_final_neoantigen.tsv"
	rnaeditor_path="software/RNAEditor/"
	time.sleep(1)
	print "Read and parse parameters done..."
	print "Check reference file path and input file path..."
	#####check input file,tool path and reference file#####
	if os.path.exists(rna_fastq_1_path) and os.path.exists(rna_fastq_2_path):
		MODE="PE"
		print "Paired-end RNA sequence provided..."
	elif os.path.exists(rna_fastq_1_path) and (not os.path.exists(rna_fastq_2_path)):
		MODE="SE"
		print "Single-end RNA sequence provided..."
	else:
		print "please check your input fastq file!"
		os._exit(1)
	if os.path.exists(opitype_fold):
		print "check opitype path done."
	else:
		print "please check your opitype path!"
		os._exit(1)
	if os.path.exists(star_path):
		print "check star path done."
	else:
		print "please check your star path!"
		os._exit(1)	
	if os.path.exists(stringtie_path):
		print "check stringtie path done."
	else:
		print "please check your stringtie path!"
		os._exit(1)	
	if os.path.exists(samtools_path):
		print "check samtools path done."
	else:
		print "please check your samtools path!"
		os._exit(1)	
	if os.path.exists(java_picard_path):
		print "check picard path done."
	else:
		print "please check your picard path!"
		os._exit(1)	
	if os.path.exists(GATK_path):
		print "check GATK path done."
	else:
		print "please check your GATK path!"
		os._exit(1)	
	if os.path.exists(vep_path):
		print "check vep path done."
	else:
		print "please check your vep path!"
		os._exit(1)	
	if os.path.exists(vep_cache):
		print "check vep cache path done."
	else:
		print "please check your vep cache path!"
		os._exit(1)		
	if os.path.exists(REFERENCE):
		print "check REFERENCE file path done."
	else:
		print "please check your REFERENCE file path!"
		os._exit(1)	
	time.sleep(1)
	#####check output directory###
	print "check output directory"
	if not os.path.exists(output_fold):
		os.mkdir(output_fold)
	if not os.path.exists(netmhc_out_fold):
		os.mkdir(netmhc_out_fold)
	if not os.path.exists(netctl_out_fold):
		os.mkdir(netctl_out_fold)
	if not os.path.exists(alignment_out_fold):
		os.mkdir(alignment_out_fold)
	if not os.path.exists(logfile_out_fold):
		os.mkdir(logfile_out_fold)
	if not os.path.exists(clean_fastq_fold):
		os.mkdir(clean_fastq_fold)	
	print "ALL file paths are correct!"	
	time.sleep(5)
	print "start fastq quality control"
	processes_0=[]
	if MODE=="PE":
		q1=multiprocessing.Process(target=read_trimmomatic_PE,args=(rna_fastq_1_path,rna_fastq_2_path,trimmomatic_path,adapter_path_PE,tumor_fastq_prefix,logfile_out_fold,CPU,))
		processes_0.append(q1)
	elif MODE=="SE":
		q1=multiprocessing.Process(target=read_trimmomatic_SE,args=(rna_fastq_1_path,trimmomatic_path,adapter_path_SE,tumor_fastq_prefix,logfile_out_fold,CPU,))
		processes_0.append(q1)		
	for p in processes_0:
		p.daemon = True
		p.start()
	for p in processes_0:
		p.join()
	print "Mapping and HLA typing!"
	processes_1=[]
	if MODE=="PE":
		h1=multiprocessing.Process(target=mapping_PE,args=(tumor_fastq_clean_1,tumor_fastq_clean_2,CPU,alignment_out_fold,prefix,star_path,star_index_path,stringtie_path,gtf_path,java_picard_path,GATK_path,REFERENCE,indels,rnaeditor_path,))
		processes_1.append(h1)
		h2=multiprocessing.Process(target=hlatyping_pe,args=(rna_fastq_1_path,rna_fastq_2_path,opitype_fold,opitype_out_fold,opitype_ext,prefix,))
		processes_1.append(h2)
	elif MODE=="SE":
		h1=multiprocessing.Process(target=mapping_SE,args=(tumor_fastq_clean_1,CPU,alignment_out_fold,prefix,star_path,star_index_path,stringtie_path,gtf_path,java_picard_path,GATK_path,REFERENCE,indels,rnaeditor_path,))
		processes_1.append(h1)
		h2=multiprocessing.Process(target=hlatyping_se,args=(rna_fastq_1_path,opitype_fold,opitype_out_fold,opitype_ext,prefix,))
		processes_1.append(h2)	
	for h in processes_1:
		h.daemon = True
		h.start()
	for h in processes_1:
		h.join()
	print "Detecting RNA editing sites!"
	hla_str=open(opitype_out_fold+'/'+prefix+"_optitype_hla_type").readlines()[0]
	print hla_str
	processes_2=[]
	j1=multiprocessing.Process(target=neo_detection,args=(alignment_out_fold,prefix,vep_path,vep_cache,human_peptide_path,netMHCpan_path,hla_str,netchop_path))
	processes_2.append(j1)
	for j in processes_2:
		j.daemon = True
		j.start()
	for j in processes_2:
		j.join()
	processes_3=[]
	t1=multiprocessing.Process(target=neo_score,args=(editing_gvf_file,neo_file,final_neo_file,cf_hy_model_9,cf_hy_model_10,cf_hy_model_11,))
	processes_3.append(t1)
	for t in processes_3:
		t.daemon = True
		t.start()
	for t in processes_3:
		t.join()
