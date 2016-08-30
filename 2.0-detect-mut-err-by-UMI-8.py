#!/usr/bin/env python
# encoding: utf-8
"""
2.0-detect-mut-err-by-UMI8.py

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


def cdhit_right_reads(right_primer_fasta):
	head, tail 	= os.path.splitext(right_primer_fasta)
	barcode 	= head.split("/")[-1].split("_")[4]
	round_index = 1
	cmd = "cd-hit-est -i %s -o %s_%s -c 0.95 -n 10 -d 0 -M 0 -T 0 -g 1 -p 1"%(right_primer_fasta, head, round_index)
	out_temp = tempfile.SpooledTemporaryFile(bufsize = 10*1000)
	fileno = out_temp.fileno()
	obj = subprocess.Popen(cmd, stdout = fileno, stderr = fileno, shell = True)
	obj.wait()
	out_temp.seek(0)
	lines = out_temp.readlines()
	print lines
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
	pool_size = multiprocessing.cpu_count()
	"""
	#'''
	# Step1: Split same UMI 8 reads to same file
	os.system("rm %s/justfy_primer_and_group_*.sh" %(prj_tree.jobs))
	os.system("rm %s/errput*" %(prj_tree.jobs))
	os.system("rm %s/output*" %(prj_tree.jobs))
	os.system("rm %s/*.fasta"%(prj_tree.clustal_fasta))
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
		for index, (key, value) in enumerate(sorted(result_dict.items())):
			#if index <= 100:
			pickle_file = '%s/%s_justfy_primer_and_group_dump_%s'%(prj_tree.tmp, prj_name, index)
			pickle_file_handle = open(pickle_file, 'wb')
			dump_tuple = (key,value)
			pickle.dump(dump_tuple, pickle_file_handle)
			pickle_file_handle.close()
			#justfy_primer_and_group(key,value,origin_fasta_dict, recombanation_info_dict, UMI_length)
			prepare_justfy_primer_and_group_pbs(prj_name, prj_tree, pickle_file)
			#group_pool.apply_async(justfy_primer_and_group, args= (key, value, shared_origin_fasta_dict, shared_recombanation_info_dict, UMI_length))
			
			
				
	pjobs = glob.glob("%s/justfy_primer_and_group_*.sh" %(prj_tree.jobs))
	pool_size = multiprocessing.cpu_count()
	pool = Pool(processes = pool_size-1)
	print "BSUB parallel jobs..."
	pjobs_ids = pool.map_async(bsub_jobs, pjobs).get(int(30*len(pjobs)))
	print "pjobs: %s pjob has been submited."%len(pjobs_ids)

	print "Waiting for all subprocesses done..."
	pool.close()
	pool.join()
	check_jobs_done(prj_name, prj_tree, "justfy_primer_and_group", pjobs_ids)
	print 'All subprocesses done.'
	#'''
	
	#'''
	#Step2: caculate composition of group
	UMI_lengths = ["_8"]
	for UMI_length in UMI_lengths:
		writer = csv.writer(open("%s/%s_cut_berfore_UMI%s_group_composition.txt"%(prj_tree.figure, prj_name, UMI_length), "w"), delimiter = "\t")
		writer.writerow(["Chain type","Barcode", "#in group","#notin group","#right primer","#wrong primer","#Total group","#Total primer"])
		in_group_files     = glob.glob("%s/%s_*_cut_berfore_UMI%s_in_group.fasta"%(prj_tree.clustal_fasta, prj_name, UMI_length))
		for index, in_group_file in enumerate(in_group_files):
			if index % 1000 == 0:
				print "Processing %s ..."%in_group_file
			
			barcode = in_group_file.split("/")[-1].split("_")[4]
			notin_group_file  = "%s/%s_%s_cut_berfore_UMI%s_notin_group.fasta"%(prj_tree.clustal_fasta, prj_name, barcode, UMI_length)
			right_primer_file = "%s/%s_%s_cut_berfore_UMI%s_right_primer.fasta"%(prj_tree.clustal_fasta, prj_name, barcode, UMI_length)
			wrong_primer_file = "%s/%s_%s_cut_berfore_UMI%s_wrong_primer.fasta"%(prj_tree.clustal_fasta, prj_name, barcode, UMI_length)
			in_group_dict = SeqIO.index(in_group_file, "fasta")
			notin_group_dict = SeqIO.index(notin_group_file, "fasta")
			right_primer_dict = SeqIO.index(right_primer_file, "fasta")
			wrong_primer_dict = SeqIO.index(wrong_primer_file, "fasta")
			try:
				chain_type = str([x for x in in_group_dict.keys()][0]).split("_")[0]
			except:
				chain_type = "No_type"
			if (len(wrong_primer_dict) + len(right_primer_dict)) != (len(in_group_dict) + len(notin_group_dict)):
				print "Warning!","primer:%s, group: %s"%((len(wrong_primer_dict) + len(right_primer_dict)),(len(in_group_dict) + len(notin_group_dict)))
				#sys.exit(0)
			#print "primer:%s, group: %s"%((len(wrong_primer_dict) + len(right_primer_dict)),(len(in_group_dict) + len(notin_group_dict)))
			writer.writerow([chain_type, barcode, len(in_group_dict), len(notin_group_dict), len(right_primer_dict), len(wrong_primer_dict)])
	#'''
	
	"""
	'''#Don't use this module, we think VJ recomb info can't sort diff Bcell which have same UMI
	#Step 4: detect err
	group_type = "right_primer"
	UMI_lengths = ["_8"]
	UMI_length = "_8"
	recomb_file = csv.reader(open("%s/%s_get_recombanation_info.txt"%(prj_tree.igblast_data, prj_name),"r"),delimiter = "\t")
	recomb_dict = {}
	for line in recomb_file:
		read_id = line[0].replace(' ', '')
		if len(line) > 1 and ("K" in line[1] or "L" in line[1]):
			recomb_dict[read_id] = (line[1], line[2])
		if len(line) > 1 and "H" in line[1] :
			recomb_dict[read_id] = (line[1], line[3])
	#right_primer_fastas = ["/zzh_gpfs02/yanmingchen/HJT-PGM/PGM_UMI_121_20160624/clustal_fasta/PGM_UMI_121_20160624_TTTTTTTT_cut_berfore_UMI_8_right_primer.fasta"]
	right_primer_fastas = glob.glob("%s/%s_*_cut_berfore_UMI%s_%s.fasta" %(prj_tree.clustal_fasta, prj_name, UMI_length,  group_type))
	for right_primer_fasta in right_primer_fastas:
		print "Processing %s " %right_primer_fasta
		umi = right_primer_fasta.split("/")[-1].split("_")[4]
		same_umi_recomb_reads_outfile = csv.writer(open("%s/%s_%s_recomb_reads.txt"%(prj_tree.clustal_fasta, prj_name, umi),"w"), delimiter = "\t")
		same_umi_recomb_reads_outfile.writerow(["Index_num", "reads_num", "V", "J", "reads_ids"])
		split_same_recomb_reads_dict = {}
		right_primer_fasta_dict = SeqIO.index(right_primer_fasta, "fasta")
		for record in right_primer_fasta_dict.values():
			record_id = record.id.split("_")[-1]
			split_same_recomb_reads_dict.setdefault(recomb_dict[record_id], []).append(record.id)
			
		for index, (VJ_recomb, same_recomb_reads_ids) in enumerate(split_same_recomb_reads_dict.items()):
			
			same_umi_recomb_reads_fasta_outfile = open("%s/%s_%s_recomb_reads_index%s_num%s.fasta"%(prj_tree.clustal_fasta, prj_name, umi, index, len(same_recomb_reads_ids)),"w")
			#result = list(VJ_reomb.extend(same_recomb_reads_ids)
			VJ_recomb_list = list(VJ_recomb)
			same_umi_recomb_reads_outfile.writerow([index, len(same_recomb_reads_ids)] + VJ_recomb_list + same_recomb_reads_ids)
			for reads_id in  same_recomb_reads_ids:
				SeqIO.write(right_primer_fasta_dict[reads_id], same_umi_recomb_reads_fasta_outfile, "fasta")
	'''
	
	
	#step3: clustal
	'''
	os.system("rm %s/clustal_*.sh" %(prj_tree.jobs))
	group_type = "right_primer"
	UMI_lengths = ["_8"]
	UMI_length = "_8"
	for UMI_length in UMI_lengths:
		prepare_clustal_jobs_normal(prj_name, prj_tree, UMI_length, group_type)
	clustal_jobs = glob.glob("%s/clustal_*.sh" %(prj_tree.jobs))
	clustal_pool = Pool()
	clustal_jobs_ids = clustal_pool.map_async(bsub_jobs, clustal_jobs).get(int(30*len(clustal_jobs)))
	print "Cluatal_jobs: %s clustal_job has been submited."%len(clustal_jobs_ids)
	
	print "Waiting for all subprocesses done..."
	clustal_pool.close()
	clustal_pool.join()
	check_jobs_done(prj_name, prj_tree, "clustal", clustal_jobs_ids)
	print 'All subprocesses done.'
	'''
	#"""

	#Step4: Consensus sequence
	
	#"""# No use vdj assign information
	group_type = "right_primer"
	UMI_lengths = ["_8"]
	UMI_length = "_8"
	clustal_alns = ["/zzh_gpfs02/yanmingchen/HJT-PGM/PGM_UMI_121_20160624/clustal_fasta/PGM_UMI_121_20160624_TTTTTTTA_cut_berfore_UMI_8_in_group.aln"]
	#clustal_alns = glob.glob("%s/%s_*_cut_berfore_UMI%s_%s.aln" %(prj_tree.clustal_fasta, prj_name, UMI_length,  group_type))
	#'''
	#Step1: Consensus sequence
	consensus_seq_fasta, consensus_seq_num = open("%s/%s_cut_berfore_UMI%s_%s_consensus.fasta"%(prj_tree.analysis, prj_name, UMI_length,  group_type), "w"), 0
	consensusX_seq_fasta, consensusX_seq_num = open("%s/%s_cut_berfore_UMI%s_%s_consensusX.fasta"%(prj_tree.analysis, prj_name, UMI_length,  group_type), "w"), 0
	consensusX_trim_seq_fasta, consensusX_trim_seq_num = open("%s/%s_cut_berfore_UMI%s_%s_consensusX_trim.fasta"%(prj_tree.analysis, prj_name, UMI_length,  group_type), "w"), 0

	for index, clustal_aln in enumerate(clustal_alns):
		barcode = clustal_aln.split("/")[-1].split("_")[4]
		if index%100 == 0 :
			print "Processing %s " %clustal_aln
		try:
			c_align = AlignIO.read(clustal_aln, 'clustal')
		except ValueError:
			continue
		sorted_reads_number = 0
		while len(c_align) - sorted_reads_number >= 3:
			remain_c_align = c_align[sorted_reads_number : ]
			summary_align = AlignInfo.SummaryInfo(remain_c_align)
			consensus_seq = summary_align.dumb_consensus(consensus_alpha = Bio.Alphabet.IUPAC.IUPACAmbiguousDNA)
			consensus_seq = str(consensus_seq)
			my_pssm = summary_align.pos_specific_score_matrix(consensus_seq)

			summary_align = AlignInfo.SummaryInfo(c_align[0:86])
			consensus_seq = summary_align.dumb_consensus(consensus_alpha = Bio.Alphabet.IUPAC.IUPACAmbiguousDNA)
			consensus_seq = str(consensus_seq)
			my_pssm = summary_align.pos_specific_score_matrix(consensus_seq)

			zs = zip(consensus_seq,my_pssm)
			nucle_percent_list = []
			for index,(ref_nucle, pssm_nums) in enumerate(zs):
				#print pssm_nums.keys()  #['A', 'C', '-', 'T', 'G']

				line = pssm_nums.values()
				line.insert(0, ref_nucle)
				pssm_num_sorted = sorted(pssm_nums.items(), key=lambda d: d[1], reverse=True)
				consensus_max_nule = pssm_num_sorted[0][0]
				consensus_seq = consensus_seq[:index] + consensus_max_nule + consensus_seq[index+1:]
				nucle_percent = pssm_nums[consensus_max_nule]/sum(pssm_nums.values()) *100
				nucle_percent_list.append(nucle_percent)
			first200_position = 0
			for index, item in enumerate(consensus_seq):
				if item != "-":
					first200_position += 1
				if first200_position == 200:
					consensus_seq = consensus_seq[:index+1]
					nucle_percent_list = nucle_percent_list[:index+1]

			nucle_percent_median = get_median_v2(nucle_percent_list)
			nucle_percent_median_index = nucle_percent_list.index(nucle_percent_median)
			percent_list_40 = nucle_percent_list[nucle_percent_median_index-20: nucle_percent_median_index+20]
			reads_number = sum(percent_list_40)/float(len(percent_list_40))
			reads_number = int(round(reads_number))
			sorted_reads_align = remain_c_align[0:reads_number]
			for index, item in enumerate(sorted_reads_align):
				print index, item
				
			sorted_reads_number += reads_number
			print reads_number
			sys.exit(0)
		while consensus_seq.startswith('X'):
			consensus_seq = consensus_seq[1:]
		while consensus_seq.endswith('X'):
			consensus_seq = consensus_seq[:-1]
		if "X" not in consensus_seq:
			consensus = SeqRecord_gernerator("%s_%s_%s"%(prj_name, barcode, VJ_recomb_index + "_" + reads_depth), str(consensus_seq), "")
			SeqIO.write(consensus, consensusX_trim_seq_fasta, "fasta")
			consensusX_trim_seq_num += 1
		else:
			consensus = SeqRecord_gernerator("%s_%s_%s"%(prj_name, barcode, VJ_recomb_index + "_" + reads_depth), str(consensus_seq), "")
			SeqIO.write(consensus, consensusX_seq_fasta, "fasta")
			consensusX_seq_num += 1
		
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
	'''
	
	
	#Step 4: detect err
	'''
	group_type = "in_group"
	UMI_lengths = ["_8"]
	UMI_length = "_8"
	clustal_alns = ["/zzh_gpfs02/yanmingchen/HJT-PGM/PGM_UMI_121_20160624/clustal_fasta/PGM_UMI_121_20160624_AAAAAAGA_cut_berfore_UMI_8_in_group.aln"]
	#clustal_alns = glob.glob("%s/%s_*_cut_berfore_UMI%s_%s.aln" %(prj_tree.clustal_fasta, prj_name, UMI_length,  group_type))
	consensus_seq_fasta = open("%s/%s_*_cut_berfore_UMI%s_%s_consensus.fasta"%(prj_tree.analysis, prj_name, UMI_length,  group_type), "w")
	for clustal_aln in clustal_alns:
		print "Processing %s " %clustal_aln
		
		
		outfile_pssm  = csv.writer(open("%s/%s_*_cut_berfore_UMI%s_%s_pssm.txt"%(prj_tree.analysis, prj_name, UMI_length,  group_type),"w"),delimiter = "\t")
		outfile_error = csv.writer(open("%s/%s_*_cut_berfore_UMI%s_%s_error.txt"%(prj_tree.analysis, prj_name, UMI_length,  group_type),"w"),delimiter = "\t")
		
		
		
		# get an alignment object from a Clustalw alignment output
		c_align = AlignIO.read(clustal_aln, 'clustal')
		summary_align = AlignInfo.SummaryInfo(c_align)
		consensus = summary_align.dumb_consensus(consensus_alpha = Bio.Alphabet.IUPAC.IUPACAmbiguousDNA)
		SeqIO.write(consensus, )
		
		my_pssm = summary_align.pos_specific_score_matrix(consensus)

		print "The consensus sequence is %s"%consensus
	
		outfile1.writerow(["Ref_nul","A","C","-","T","G"])
		zs = zip(consensus,my_pssm)
		for index,(i,j) in enumerate(zs):
			line = j.values()
			line.insert(0,i)
			outfile1.writerow(line)
	
		replace_info = summary_align.replacement_dictionary()

		outfile2.writerow(['Position','Type','Original','Later'])
		zs = zip(c_align[0].seq,c_align[1].seq)
		for index,(i,j) in enumerate(zs):
			if i != j :
				if i == '-':
					outfile2.writerow([index,'Insertion',i,j])
				elif j == '-':
					outfile2.writerow([index,'Deletion',i,j])
				else:
					outfile2.writerow([index,'Mutation',i,j])
	'''
	#"""
	