#!/usr/bin/env python
# encoding: utf-8
"""
2.2-trim-productive-reads.py

Created by Mingchen on 2015-01-19.
Copyright (c) 2015 __MyCompanyName__. All rights reserved.

Usage:


"""
import sys, os, csv, re, glob, copy, subprocess, time, multiprocessing
import traceback, tempfile
import Bio.Alphabet 
from Bio.Align import AlignInfo
from Bio import AlignIO
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


if __name__=='__main__':
	print 'Parent process %s'%os.getpid()
	prj_folder = os.getcwd()
	
	prj_tree = ProjectFolders(prj_folder)
	
	prj_name = fullpath2last_folder(prj_tree.home)
	pool_size = multiprocessing.cpu_count()
	#group_type = "right_primer"
	group_type = "in_group"
	UMI_lengths = ["_8"]
	UMI_length = "_8"

	"""
	# Trim productive reads
	origin_fasta = "%s/%s.fasta"%(prj_tree.origin, prj_name)
	assign_file = csv.reader(open("%s/%s_get_assignment_info.txt"%(prj_tree.igblast_data, prj_name),"r"),delimiter = "\t")
	recomb_file = csv.reader(open("%s/%s_get_recombanation_info.txt"%(prj_tree.igblast_data, prj_name),"r"),delimiter = "\t")
	
	productive_reads_id_list, assigned_reads_id_list = [], []
	for line in recomb_file:
		if len(line) < 2:
			continue
		if line[-2] == "Yes" or line[-2] == "No":
			if line[-1] == "+":
				read_id = line[0].replace(' ', '')
			elif line[-1] == "-":
				read_id = "reversed|"+line[0].replace(' ', '')
			assigned_reads_id_list.append(read_id)
		if line[-2] == "Yes":
			if line[-1] == "+":
				read_id = line[0].replace(' ', '')
			elif line[-1] == "-":
				read_id = "reversed|"+line[0].replace(' ', '')
			productive_reads_id_list.append(read_id)
	print len(productive_reads_id_list)
	outfile = open("%s/%s_trimed_productive_variable_region.txt"%(prj_tree.analysis, prj_name),"w")
	#outfile = open("%s/%s_trimed_variable_region.txt"%(prj_tree.analysis, prj_name),"w")
	writer = csv.writer(outfile,delimiter = "\t")
	trimed_productive_fasta = open("%s/%s_trimed_productive_variable_region.fasta"%(prj_tree.analysis, prj_name),"w")
	#trimed_productive_fasta = open("%s/%s_trimed_variable_region.fasta"%(prj_tree.analysis, prj_name),"w")
	result = ['Query_ID','query_seq','query_length','trimmed_reads','tr_length']
	origin_fasta_dict = load_fastas_in_list_v2(origin_fasta, productive_reads_id_list)
	ID, count, count1, start, end = 'start_id', 0, 0, int(), int()
	for line in assign_file:
		query_ID = line[1]
		nr_query_ID = query_ID.replace('reversed|','')
		if query_ID in productive_reads_id_list :
			if query_ID != ID:
				count1 += 1
				if (count1) % 10000 == 0:
					print "%d reads processed" %(count1)
				if 'reversed' in query_ID:
					count += 1
					a = origin_fasta_dict.get(nr_query_ID)
					rc = a.reverse_complement(id = query_ID )
					query_seq = rc.seq
				else:
					query_seq = origin_fasta_dict.get(nr_query_ID).seq
				writer.writerow(result)
				trimed_productive_read = SeqRecord_gernerator(result[0], result[3], "")
				SeqIO.write(trimed_productive_read, trimed_productive_fasta, "fasta")
				ID = query_ID
				start = int(line[8])
				#print line[10],line
				trans_start = int(line[10])
				end = int(line[9])
			else:
				end = int(line[9])
			result[0] = query_ID
			result[1] = query_seq
			result[2] = len(query_seq)
			if start-trans_start < 0:
				
				if (trans_start-1) % 3 == 0:
					result[3] = query_seq[start -1 :end]
				else:
					#Wrong! and beside correceted! result[3] = query_seq[(start - 1 + 3 - (trans_start-1)%3 -1):end]
					result[3] = query_seq[(start - 1 + 2 - ((trans_start-1)%3 -1)):end]
			else:
				if (trans_start-1) % 3 == 0:
					result[3] = query_seq[start -1 :end]
				else:
					#Wrong! and beside correceted! result[3] = query_seq[(start - 1 + 3 - (trans_start-1)%3 -1):end]
					result[3] = query_seq[(start - 1 + 2 - ((trans_start-1)%3 -1)):end]
				#print "Maybe a wrong reads which not start at first nucle.%s"%query_ID
				#result[3] = query_seq[start-(trans_start-1)-1:end]
				#result[3] = query_seq[(start - 1 + 3 - (trans_start-1)%3 -1):end]
			result[4] = len(result[3])
	writer.writerow(result)
	trimed_productive_read = SeqRecord_gernerator(result[0], result[3], "")
	SeqIO.write(trimed_productive_read, trimed_productive_fasta, "fasta")
	print "There are %s productive productive!"%(len(productive_reads_id_list))
	print 'Find %d seq and %d reversed seq...for %s'%(count1,count,fname)
	"""
	
	# Trim all assgned reads
	origin_fasta = "%s/%s.fasta"%(prj_tree.origin, prj_name)
	assign_file = csv.reader(open("%s/%s_get_assignment_info.txt"%(prj_tree.igblast_data, prj_name),"r"),delimiter = "\t")
	recomb_file = csv.reader(open("%s/%s_get_recombanation_info.txt"%(prj_tree.igblast_data, prj_name),"r"),delimiter = "\t")
	
	productive_reads_id_list, assigned_reads_id_list = [], []
	for line in recomb_file:
		if len(line) < 2:
			continue
		if line[-2] == "Yes" or line[-2] == "No":
			if line[-1] == "+":
				read_id = line[0].replace(' ', '')
			elif line[-1] == "-":
				read_id = "reversed|"+line[0].replace(' ', '')
			assigned_reads_id_list.append(read_id)
	print len(assigned_reads_id_list)
	outfile = open("%s/%s_trimed_variable_region.txt"%(prj_tree.analysis, prj_name),"w")
	writer = csv.writer(outfile,delimiter = "\t")
	trimed_assigned_fasta = open("%s/%s_trimed_variable_region.fasta"%(prj_tree.analysis, prj_name),"w")
	result = ['Query_ID','query_seq','query_length','trimmed_reads','tr_length']
	origin_fasta_dict = load_fastas_in_list_v2(origin_fasta, assigned_reads_id_list)
	ID, count, count1, start, end = 'start_id', 0, 0, int(), int()
	for line in assign_file:
		query_ID = line[1]
		nr_query_ID = query_ID.replace('reversed|','')
		if query_ID in assigned_reads_id_list :
			if query_ID != ID:
				count1 += 1
				if (count1) % 10000 == 0:
					print "%d reads processed" %(count1)
				if 'reversed' in query_ID:
					count += 1
					a = origin_fasta_dict.get(nr_query_ID)
					rc = a.reverse_complement(id = query_ID )
					query_seq = rc.seq
				else:
					query_seq = origin_fasta_dict.get(nr_query_ID).seq
				writer.writerow(result)
				trimed_assigned_read = SeqRecord_gernerator(result[0], result[3], "")
				SeqIO.write(trimed_assigned_read, trimed_assigned_fasta, "fasta")
				ID = query_ID
				start = int(line[8])
				#print line[10],line
				trans_start = int(line[10])
				end = int(line[9])
			else:
				end = int(line[9])
			result[0] = query_ID
			result[1] = query_seq
			result[2] = len(query_seq)
			if start-trans_start < 0:
				
				if (trans_start-1) % 3 == 0:
					result[3] = query_seq[start -1 :end]
				else:
					#Wrong! and beside correceted! result[3] = query_seq[(start - 1 + 3 - (trans_start-1)%3 -1):end]
					result[3] = query_seq[(start - 1 + 2 - ((trans_start-1)%3 -1)):end]
			else:
				if (trans_start-1) % 3 == 0:
					result[3] = query_seq[start -1 :end]
				else:
					#Wrong! and beside correceted! result[3] = query_seq[(start - 1 + 3 - (trans_start-1)%3 -1):end]
					result[3] = query_seq[(start - 1 + 2 - ((trans_start-1)%3 -1)):end]
				#print "Maybe a wrong reads which not start at first nucle.%s"%query_ID
				#result[3] = query_seq[start-(trans_start-1)-1:end]
				#result[3] = query_seq[(start - 1 + 3 - (trans_start-1)%3 -1):end]
			result[4] = len(result[3])
	writer.writerow(result)
	trimed_assigned_read = SeqRecord_gernerator(result[0], result[3], "")
	SeqIO.write(trimed_assigned_read, trimed_assigned_fasta, "fasta")
	print "There are %s assigned reads!"%(len(assigned_reads_id_list))
	print 'Find %d seq and %d reversed seq...for %s'%(count1,count,fname)	
				