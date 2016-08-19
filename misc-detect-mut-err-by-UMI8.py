#!/usr/bin/env python
# encoding: utf-8
"""
misc-detect-mut-err-by-UMI8.py

Created by Mingchen on 2015-01-19.
Copyright (c) 2015 __MyCompanyName__. All rights reserved.

Usage: misc-get-pattern-loc.py -i infile -r referencefile -o outfile


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
	#clustal_alns = ["/zzh_gpfs02/yanmingchen/HJT-PGM/PGM_UMI_121_20160624/clustal_fasta/PGM_UMI_121_20160624_TTTTTTTA_cut_berfore_UMI_8_in_group.aln"]
	clustal_alns = glob.glob("%s/%s_*_cut_berfore_UMI%s_%s.aln" %(prj_tree.clustal_fasta, prj_name, UMI_length,  group_type))
	'''
	#Step1: Consensus sequence
	consensus_seq_fasta, consensus_seq_num = open("%s/%s_cut_berfore_UMI%s_%s_consensus.fasta"%(prj_tree.analysis, prj_name, UMI_length,  group_type), "w"), 0
	consensusX_seq_fasta, consensusX_seq_num = open("%s/%s_cut_berfore_UMI%s_%s_consensusX.fasta"%(prj_tree.analysis, prj_name, UMI_length,  group_type), "w"), 0
	for index, clustal_aln in enumerate(clustal_alns):
		barcode = clustal_aln.split("/")[-1].split("_")[4]
		if index%100 == 0 :
			print "Processing %s " %clustal_aln
		try:
			c_align = AlignIO.read(clustal_aln, 'clustal')
		except ValueError:
			continue
		summary_align = AlignInfo.SummaryInfo(c_align)
		consensus_seq = summary_align.dumb_consensus(consensus_alpha = Bio.Alphabet.IUPAC.IUPACAmbiguousDNA)
		consensus = SeqRecord_gernerator("%s_%s_%s"%(prj_name, barcode, index), str(consensus_seq), "")
		
		#print consensus_seq
		if "X" not in consensus_seq:
			SeqIO.write(consensus, consensus_seq_fasta, "fasta")
			consensus_seq_num += 1
		else:
			SeqIO.write(consensus, consensusX_seq_fasta, "fasta")
			consensusX_seq_num += 1
	print consensus_seq_num, consensusX_seq_num, len(clustal_alns)
	'''
	#Step2: Error
	for index, clustal_aln in enumerate(clustal_alns):
		barcode = clustal_aln.split("/")[-1].split("_")[4]
		if index%100 == 0 :
			print "Processing %s " %clustal_aln
		try:
			c_align = AlignIO.read(clustal_aln, 'clustal')
		except ValueError:
			continue
		summary_align = AlignInfo.SummaryInfo(c_align)
		consensus_seq = summary_align.dumb_consensus(consensus_alpha = Bio.Alphabet.IUPAC.IUPACAmbiguousDNA)
		consensus = SeqRecord_gernerator("%s_%s_%s"%(prj_name, barcode, index), str(consensus_seq), "")
	
		my_pssm = summary_align.pos_specific_score_matrix(consensus)
		outfile_pssm  = csv.writer(open("%s/%s_*_cut_berfore_UMI%s_%s_pssm.txt"%(prj_tree.analysis, prj_name, UMI_length,  group_type),"w"),delimiter = "\t")
		outfile_error = csv.writer(open("%s/%s_*_cut_berfore_UMI%s_%s_error.txt"%(prj_tree.analysis, prj_name, UMI_length,  group_type),"w"),delimiter = "\t")
		outfile_error.writerow(['Position','Type','Original','Later','Percent'])
		#for line in my_pssm:
			#print line.keys(), line.values()
		if "X" in consensus_seq:
			print my_pssm
			sys.exit(0)
		zs = zip(consensus,my_pssm)
		for index,(ref_nucle, pssm_nums) in enumerate(zs):
			line = pssm_nums.values()
			line.insert(0, ref_nucle)
			if line[1] != 0 and ref_nucle != "A":
				mut_list.append([index,'Mutation',ref_nucle,'A'])
			
			if line[2] != 0 and ref_nucle != "C":
				mut_list.append([index,'Mutation',ref_nucle,'C'])
			
			if line[4] != 0 and ref_nucle != "T":
				mut_list.append([index,'Mutation',ref_nucle,'T'])
			
			if line[5] != 0 and ref_nucle != "G":
				mut_list.append([index,'Mutation',ref_nucle,'G'])
			if line[3] != 0:
				del_list.append([index,'Deletion',ref_nucle,'_'])
				pos_list.append(index)
			
			outfile_error.writerow([index,'Type','Original','Later'])
		for read_index in range(0, len(c_align)):
			zs = zip(str(consensus_seq),c_align[read_index].seq)
			print len(str(consensus_seq)), len(c_align[read_index].seq)
			print zs
		
			for index,(i,j) in enumerate(zs):
				if i != j :
					if i == '-':
						outfile2.writerow([index,'Insertion',i,j])
					elif j == '-':
						outfile2.writerow([index,'Deletion',i,j])
					else:
						outfile2.writerow([index,'Mutation',i,j])
		#if len(c_align) >= 3 and len(c_align) <= 10:
			 