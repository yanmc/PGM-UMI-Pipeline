#!/usr/bin/env python
# encoding: utf-8
"""
2.0-detect-mut-err-by-UMI8.py

Created by Mingchen on 2015-01-19.
Copyright (c) 2015 __MyCompanyName__. All rights reserved.

Usage: misc-get-pattern-loc.py -i infile -r referencefile -o outfile


"""
import sys, os, csv, re, glob, copy, subprocess, time
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from multiprocessing import Pool
from bsub import bsub
from mytools import *

def load_recombanation_info_dict(handle):
	recombanation_info_file = csv.reader(open(handle, "rU"), delimiter="\t")
	recombanation_info_dict = {}
	for line in recombanation_info_file:
		line[0] = line[0].replace(" ", "")
		if len(line) > 1 and "H" in line[1]:
			recombanation_info_dict[line[0]] = "H"
		
		elif len(line) > 1 and "K" in line[1]:
			recombanation_info_dict[line[0]] = "K"
		
		elif len(line) > 1 and "L" in line[1]:
			recombanation_info_dict[line[0]] = "L"
		else:
			recombanation_info_dict[line[0]] = 0
	
	return recombanation_info_dict

def write_same_UMI_file(result_dict, infile):
	infile_index = retrieve_name_body(infile).split("-")[0].split("_")[-1]
	handle = csv.reader(open(infile, "rU"), delimiter='\t')
	for line in handle:
		UMI = line[0]
		reads_id = line[1:]
		result_dict.setdefault(UMI, []).extend(reads_id)
	return result_dict

def bsub_jobs(job):
	cmd = '%s <  %s'%('bsub', job)
	head, tail 	= os.path.splitext(job)	
	p=subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)  
	while True:
		buff = p.stdout.readline()
		if ">" in buff and "<" in buff:
	 		igblast_job_id = buff.split(">")[0].split("<")[-1]
		if buff == '' and p.poll() != None:
			break
	return igblast_job_id

def check_jobs_done(app, igblast_job_ids):
	log_file = "%s/%s_pbs.log"%(prj_tree.logs, app)
	log_file_handle = open(log_file, "w")
	for igblast_job_id in igblast_job_ids:
		while os.path.exists("%s/output_%s"%(prj_tree.home, igblast_job_id)) == False:
			log_file_handle.write("Waiting for job_%s...%s\n"%(igblast_job_id,time.ctime()))
			time.sleep(1)
		
		errput, output = "%s/errput_%s"%(prj_tree.home, igblast_job_id), "%s/output_%s"%(prj_tree.home, igblast_job_id)
		IgBLAST_log = True
		while IgBLAST_log:
			output_log = open(output, "rU")
			for line in output_log.readlines():
				if line.replace('\n','') == "Successfully completed.":
					log_file_handle.write("job_%s Successfully completed.\n"%igblast_job_id)
					IgBLAST_log = False
					os.system("mv %s %s"%(errput, prj_tree.logs))
					os.system("mv %s %s"%(output, prj_tree.logs))
			output_log.close()
if __name__=='__main__':
	print 'Parent process %s'%os.getpid()
	prj_folder = os.getcwd()
	
	prj_tree = ProjectFolders(prj_folder)
	
	prj_name = fullpath2last_folder(prj_tree.home)
	'''
	# Step1: Split same UMI 8 reads to same file
	print "Loading origin fasta dict......"
	
	os.system("cat %s/%s*_cut_berfore_UMI_*.fasta >  %s/%s_cut_berfore_UMI.fasta"%(prj_tree.split, prj_name, prj_tree.origin, prj_name))
	origin_fasta_dict = SeqIO.index("%s/%s_cut_berfore_UMI.fasta"%(prj_tree.origin, prj_name), "fasta")

	print "reads number:%s"%(len(origin_fasta_dict))
	
	recombanation_info_file = "%s/%s_get_recombanation_info.txt"%(prj_tree.igblast_data, prj_name)
	recombanation_info_dict = load_recombanation_info_dict(recombanation_info_file)
	
	UMI_lengths = ["","_6_8","_8"]
	for UMI_length in UMI_lengths:
		infiles = glob.glob("%s/%s*-get-primer-reads_UMI%s.txt"%(prj_tree.data, prj_name, UMI_length))
		result_dict = {}
		for infile in infiles:
			print "Processing %s ..."%infile
			result_dict = write_same_UMI_file(result_dict, infile)
		total_reads_number = 0
		for index, (key, value) in enumerate(result_dict.items()):
			print "Write No.%s UMI, it has %s reads."%(index, len(value))
			total_reads_number += len(value)
			writer = open("%s/%s_%s_cut_berfore_UMI%s_in_group.fasta"%(prj_tree.clustal_fasta, prj_name, key, UMI_length), "w")
			writer2 = open("%s/%s_%s_cut_berfore_UMI%s_notin_group.fasta"%(prj_tree.clustal_fasta, prj_name, key, UMI_length), "w")
			writer3 = open("%s/%s_%s_cut_berfore_UMI%s_right_primer.fasta"%(prj_tree.clustal_fasta, prj_name, key, UMI_length), "w")
			writer4 = open("%s/%s_%s_cut_berfore_UMI%s_wrong_primer.fasta"%(prj_tree.clustal_fasta, prj_name, key, UMI_length), "w")
			recombanation_info_list = []
			for ID in value:
				recombanation_info_dict.setdefault(ID, 0)
				try :
					recombanation_info_list.append(recombanation_info_dict[ID])
				except:
					pass
			UMI_group_type = max(map(lambda x: (recombanation_info_list.count(x), x), recombanation_info_list))[1]
			#print "UMI_group_type: %s"%UMI_group_type
			for ID in value:
				#justify primer right or wrong
				if recombanation_info_dict[ID] == "H":
					try:
						SeqIO.write(origin_fasta_dict["IGG_"+ID], writer3, "fasta")
					except:
						pass
					try:
						SeqIO.write(origin_fasta_dict["IGM_"+ID], writer3, "fasta")
					except:
						pass
				
				elif recombanation_info_dict[ID] == "K":
					try:
						SeqIO.write(origin_fasta_dict["IGK_"+ID], writer3, "fasta")
					except:
						pass
				elif recombanation_info_dict[ID] == "L":
					try:
						SeqIO.write(origin_fasta_dict["IGL_"+ID], writer3, "fasta")
					except:
						pass
				else:
					try:
						SeqIO.write(origin_fasta_dict["IGG_"+ID], writer4, "fasta")
					except:
						pass
					try:
						SeqIO.write(origin_fasta_dict["IGM_"+ID], writer4, "fasta")
					except:
						pass
					try:
						SeqIO.write(origin_fasta_dict["IGK_"+ID], writer4, "fasta")
					except:
						pass
					try:
						SeqIO.write(origin_fasta_dict["IGL_"+ID], writer4, "fasta")
					except:
						pass
				#justify in group or not
				if recombanation_info_dict[ID] == UMI_group_type == "H":
					try:
						SeqIO.write(origin_fasta_dict["IGG_"+ID], writer, "fasta")
					except:
						pass
					try:
						SeqIO.write(origin_fasta_dict["IGM_"+ID], writer, "fasta")
					except:
						pass
				elif recombanation_info_dict[ID] == UMI_group_type == "K":
					try:
						SeqIO.write(origin_fasta_dict["IGK_"+ID], writer, "fasta")
					except:
						pass
				elif recombanation_info_dict[ID] == UMI_group_type == "L":
					try:
						SeqIO.write(origin_fasta_dict["IGL_"+ID], writer, "fasta")
					except:
						pass
				elif recombanation_info_dict[ID] == UMI_group_type == 0:
					try:
						SeqIO.write(origin_fasta_dict["IGG_"+ID], writer, "fasta")
					except:
						pass
					try:
						SeqIO.write(origin_fasta_dict["IGM_"+ID], writer, "fasta")
					except:
						pass
					try:
						SeqIO.write(origin_fasta_dict["IGK_"+ID], writer, "fasta")
					except:
						pass
					try:
						SeqIO.write(origin_fasta_dict["IGL_"+ID], writer, "fasta")
					except:
						pass
				else:
					#print "%s is an error in UMI group %s !"%(ID, key)
					try:
						SeqIO.write(origin_fasta_dict["IGG_"+ID], writer2, "fasta")
					except:
						pass
					try:
						SeqIO.write(origin_fasta_dict["IGM_"+ID], writer2, "fasta")
					except:
						pass
					try:
						SeqIO.write(origin_fasta_dict["IGK_"+ID], writer2, "fasta")
					except:
						pass
					try:
						SeqIO.write(origin_fasta_dict["IGL_"+ID], writer2, "fasta")
					except:
						pass
		print 'Split UMI to same file: done, There are %s reads been processed, Raw data has %s reads.'%(total_reads_number, len(origin_fasta_dict))
	'''
	#step2: clustal
	"""
	#UMI_lengths = ["","_6_8","_8"]
	UMI_lengths = ["_8"]
	for UMI_length in UMI_lengths:
		prepare_clustal_jobs_normal(prj_name, prj_tree, UMI_length)
	clustal_jobs = glob.glob("%s/clustal_*.sh" %(prj_tree.jobs))
	clustal_pool = Pool()
	clustal_jobs_ids = clustal_pool.map_async(bsub_jobs, clustal_jobs).get(int(30*len(clustal_jobs)))
	print "Cluatal_jobs: %s clustal_job has been submited."%len(clustal_jobs_ids)
	
	print "Waiting for all subprocesses done..."
	clustal_pool.close()
	clustal_pool.join()
	check_jobs_done("clustal", clustal_jobs_ids)
	print 'All subprocesses done.'
	"""
	#Step3: caculate composition of group
	UMI_lengths = ["_8"]
	for UMI_length in UMI_lengths:
		writer = csv.writer(open("%s/%s_cut_berfore_UMI%s_group_composition.txt"%(prj_tree.figure, prj_name, UMI_length), "w"), delimiter = "\t")
		writer.writerow(["Chain type","Barcode", "#in group","#notin group","#right primer","#wrong primer","#Total group","#Total primer"])
		in_group_files     = glob.glob("%s/%s_*_cut_berfore_UMI%s_in_group.fasta"%(prj_tree.clustal_fasta, prj_name, UMI_length))
		for in_group_file in in_group_files:
			print "Processing %s ..."%in_group_file
			barcode = in_group_file.split("/")[-1].split("_")[4]
			notin_group_file  = "%s/%s_%s_cut_berfore_UMI%s_notin_group.fasta"%(prj_tree.clustal_fasta, prj_name, barcode, UMI_length)
			right_primer_file = "%s/%s_%s_cut_berfore_UMI%s_right_primer.fasta"%(prj_tree.clustal_fasta, prj_name, barcode, UMI_length)
			wrong_primer_file = "%s/%s_%s_cut_berfore_UMI%s_wrong_primer.fasta"%(prj_tree.clustal_fasta, prj_name, barcode, UMI_length)
			in_group_dict = SeqIO.index(in_group_file, "fasta")
			notin_group_dict = SeqIO.index(notin_group_file, "fasta")
			right_primer_dict = SeqIO.index(right_primer_file, "fasta")
			wrong_primer_dict = SeqIO.index(wrong_primer_file, "fasta")
			chain_type = str([in_group_dict.keys()][0]).split("_")[0]
			if (len(wrong_primer_dict) + len(right_primer_dict)) != (len(in_group_dict) + len(notin_group_dict)):
				print "Warning!","primer:%s, group: %s"%((len(wrong_primer_dict) + len(right_primer_dict)),(len(in_group_dict) + len(notin_group_dict)))
				sys.exit(0)
			writer.writerow([chain_type, barcode, len(in_group_dict), len(notin_group_dict), len(right_primer_dict), len(wrong_primer_dict)])
			
			