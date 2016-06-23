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

if __name__=='__main__':
	print 'Parent process %s'%os.getpid()
	prj_folder = os.getcwd()
	
	prj_tree = ProjectFolders(prj_folder)
	
	prj_name = fullpath2last_folder(prj_tree.home)
	
	# Split same UMI 8 reads to same file
	print "Loading origin fasta dict......"
	
	os.system("cat %s/%s*_cut_berfore_UMI_*.fasta >  %s/%s_cut_berfore_UMI.fasta"%(prj_tree.split, prj_name, prj_tree.origin, prj_name))
	origin_fasta_dict = SeqIO.index("%s/%s_cut_berfore_UMI.fasta"%(prj_tree.origin, prj_name), "fasta")
	print "reads number:%s"%(len(origin_fasta_dict))
	
	recombanation_info_file = "%s/%s_get_recombanation_info.txt"%(prj_tree.igblast_data, prj_name)
	recombanation_info_dict = load_recombanation_info_dict(recombanation_info_file)
	
	
	infiles = glob.glob("%s/%s*-get-primer-reads_UMI_8.txt"%(prj_tree.data, prj_name))
	result_dict = {}
	for infile in infiles:
		print "Processing %s ..."%infile
		result_dict = write_same_UMI_file(result_dict, infile)
	total_reads_number = 0
	for index, (key, value) in enumerate(result_dict.items()):
		print "Write No.%s UMI, it has %s reads."%(index, len(value))
		total_reads_number += len(value)
		writer = open("%s/%s_%s_cut_berfore_UMI.fasta"%(prj_tree.clustal_fasta, prj_name, key), "w")
		writer2 = open("%s/%s_%s_cut_berfore_UMI_wrong_primer.fasta"%(prj_tree.clustal_fasta, prj_name, key), "w")
		recombanation_info_list = []
		for ID in value:
			recombanation_info_dict.setdefault(ID, 0)
			try :
				recombanation_info_list.append(recombanation_info_dict[ID])
			except:
				pass
		UMI_group_type = max(map(lambda x: (recombanation_info_list.count(x), x), recombanation_info_list))[1]
		print "UMI_group_type: %s"%UMI_group_type
		for ID in value:
			
			if recombanation_info_dict[ID] == UMI_group_type == "H":
				try:
					SeqIO.write(origin_fasta_dict["IGG_"+ID], writer, "fasta")
				except:
					print "IGG type %s not in origin_fasta_dict!"%ID
				try:
					SeqIO.write(origin_fasta_dict["IGM_"+ID], writer, "fasta")
				except:
					print "IGM type %s not in origin_fasta_dict!"%ID
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
					print "IGG type %s not in origin_fasta_dict!"%ID
				try:
					SeqIO.write(origin_fasta_dict["IGM_"+ID], writer, "fasta")
				except:
					print "IGM type %s not in origin_fasta_dict!"%ID
				try:
					SeqIO.write(origin_fasta_dict["IGK_"+ID], writer, "fasta")
				except:
					print "IGK type %s not in origin_fasta_dict!"%ID
				try:
					SeqIO.write(origin_fasta_dict["IGL_"+ID], writer, "fasta")
				except:
					print "IGL type %s not in origin_fasta_dict!"%ID
			else:
				print "%s is an error in UMI group %s !"%(ID, key)
				try:
					SeqIO.write(origin_fasta_dict["IGG_"+ID], writer2, "fasta")
				except:
					print "IGG type %s not in origin_fasta_dict!"%ID
				try:
					SeqIO.write(origin_fasta_dict["IGM_"+ID], writer2, "fasta")
				except:
					print "IGM type %s not in origin_fasta_dict!"%ID
				try:
					SeqIO.write(origin_fasta_dict["IGK_"+ID], writer2, "fasta")
				except:
					print "IGK type %s not in origin_fasta_dict!"%ID
				try:
					SeqIO.write(origin_fasta_dict["IGL_"+ID], writer2, "fasta")
				except:
					print "IGL type %s not in origin_fasta_dict!"%ID
	print 'Split UMI to same file: done, There are %s reads been processed, Raw data has %s reads.'%(total_reads_number, len(origin_fasta_dict))
