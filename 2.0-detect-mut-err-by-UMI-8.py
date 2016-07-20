#!/usr/bin/env python
# encoding: utf-8
"""
2.0-detect-mut-err-by-UMI8.py

Created by Mingchen on 2015-01-19.
Copyright (c) 2015 __MyCompanyName__. All rights reserved.

Usage: misc-get-pattern-loc.py -i infile -r referencefile -o outfile


"""
import sys, os, csv, re, glob, copy, subprocess, time, multiprocessing
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from multiprocessing import Pool, Process, Manager
from bsub import bsub
from mytools import *
from misc_prepare_pbs import *
try:
    import cPickle as pickle
except ImportError:
    import pickle

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
	
if __name__=='__main__':
	print 'Parent process %s'%os.getpid()
	prj_folder = os.getcwd()
	
	prj_tree = ProjectFolders(prj_folder)
	
	prj_name = fullpath2last_folder(prj_tree.home)
	#'''
	# Step1: Split same UMI 8 reads to same file
	
	#UMI_lengths = ["","_6_8","_8"]
	UMI_lengths = ["_8"]
	for UMI_length in UMI_lengths:
		infiles = glob.glob("%s/%s*-get-primer-reads_UMI%s.txt"%(prj_tree.data, prj_name, UMI_length))
		result_dict = {}
		for infile in infiles:
			print "Processing %s ..."%infile
			result_dict = write_same_UMI_file(result_dict, infile)
		pool_size = multiprocessing.cpu_count()
		group_pool = Pool(processes = pool_size-1)
		print "Prepare parallel jobs..."
		for index, (key, value) in enumerate(result_dict.items()):
			pickle_file = '%s/%s_justfy_primer_and_group_dump_%s'%(prj_tree.tmp, prj_name, index)
			pickle_file_handle = open(pickle_file, 'wb')
			pickle.dump((key, value), pickle_file_handle)
			pickle_file_handle.close()
			#justfy_primer_and_group(key,value,origin_fasta_dict, recombanation_info_dict, UMI_length)
			prepare_justfy_primer_and_group_pbs(prj_name, prj_tree, pickle_file)
			#group_pool.apply_async(justfy_primer_and_group, args= (key, value, shared_origin_fasta_dict, shared_recombanation_info_dict, UMI_length))
		
	
	pjobs = glob.glob("%s/justfy_primer_and_group_*.sh" %(prj_tree.jobs))
	pool = Pool()
	pjobs_ids = pool.map_async(bsub_jobs, pjobs).get(int(30*len(pjobs)))
	print "pjobs: %s pjob has been submited."%len(pjobs_ids)

	print "Waiting for all subprocesses done..."
	pool.close()
	pool.join()
	check_jobs_done(prj_name, prj_tree, "justfy_primer_and_group", pjobs_ids)
	print 'All subprocesses done.'
	#'''
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
	check_jobs_done(prj_name, prj_tree, "clustal", clustal_jobs_ids)
	print 'All subprocesses done.'
	"""
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
	"""	
			