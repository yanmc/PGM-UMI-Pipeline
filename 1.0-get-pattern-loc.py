#!/usr/bin/env python
# encoding: utf-8
"""
misc-get-pattern-loc.py

Created by Mingchen on 2015-01-19.
Copyright (c) 2015 __MyCompanyName__. All rights reserved.

Usage: misc-get-pattern-loc.py -i infile -r referencefile -o outfile


"""
import sys, os, csv, re, glob, copy, subprocess, time, multiprocessing
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from multiprocessing import Pool, Process, Manager, Semaphore
import matplotlib.pyplot as plt
from mytools import *
def load_recombanation_info_dict(handle):
	recombanation_info_file = csv.reader(open(handle, "rU"), delimiter="\t")
	recombanation_info_dict = {}
	for line in recombanation_info_file:
		line[0] = line[0].replace(" ", "")
		if len(line) > 1 and "H" in line[1]:
			recombanation_info_dict[line[0]] = "H"
		
		if len(line) > 1 and "K" in line[1]:
			recombanation_info_dict[line[0]] = "K"
		
		if len(line) > 1 and "L" in line[1]:
			recombanation_info_dict[line[0]] = "L"
	return recombanation_info_dict
def get_primer_and_cut_before_UMI(handle,recombanation_info_dict):
	file_handle = open("%s"%handle, 'rU')
	fname = retrieve_name_body(handle)
	some_error_reads_file = open("%s/%s_error_reads_file.txt"%(prj_tree.analysis, fname), "w")
	cut_berfore_UMI_fasta = open("%s/%s_cut_before_UMI.fasta"%(prj_tree.tmp, fname), "w")
	outfile  = csv.writer(open("%s/%s-get-primer"%(prj_tree.data,fname),'w'),delimiter = '\t')
	outfile1 = csv.writer(open("%s/%s-get-IgG-primer"%(prj_tree.tmp,fname),'w'),delimiter = '\t')
	outfile2 = csv.writer(open("%s/%s-get-IgK-primer"%(prj_tree.tmp,fname),'w'),delimiter = '\t')
	outfile3 = csv.writer(open("%s/%s-get-IgM-primer"%(prj_tree.tmp,fname),'w'),delimiter = '\t')
	outfile4 = csv.writer(open("%s/%s-get-IgL-primer"%(prj_tree.tmp,fname),'w'),delimiter = '\t')
	#outfile1.writerow(['id','primer','start','end','UMI',"Seq","Read_type"])
	#outfile2.writerow(['id','primer','start','end','UMI',"Seq","Read_type"])
	#outfile3.writerow(['id','primer','start','end','UMI',"Seq","Read_type"])
	#outfile4.writerow(['id','primer','start','end','UMI',"Seq","Read_type"])
	IGG_primer = 'AAGACCGATGGGCCCTTGGTGGA'
	IGK_primer = 'AAGACAGATGGTGCAGCCACAGTTC'
	IGM_primer = 'AAGGGTTGGGGCGGATGCACTCC'
	IGL_primer = 'CCTTGTTGGCTTG(A|G)AGCTCCTCAGAGGAGG'
	#print "Searching...%s"%fname
	the_handle = SeqIO.parse(file_handle, 'fasta')
	for record in the_handle:
		#print "Searching...%s"%record.id
		loc1 = re.finditer('AAGACCGATGGGCCCTTGGTGGA',record.seq.tostring(),re.I)
		loc2 = re.finditer('AAGACAGATGGTGCAGCCACAGTTC',record.seq.tostring(),re.I)
		loc3 = re.finditer('AAGGGTTGGGGCGGATGCACTCC',record.seq.tostring(),re.I)
		loc4 = re.finditer('CCTTGTTGGCTTG(A|G)AGCTCCTCAGAGGAGG',record.seq.tostring(),re.I)
		try:
			Read_type = recombanation_info_dict[str(record.id)]
		except:
			Read_type = 'No_annotation'
			#print "Reads: %s is no type!"%record.id		
		for index, i in enumerate(loc1):
			if index == 0 and i.span()[0]-8 >= 0:
				#print type(i.group()),i.span()
				outfile1.writerow([record.id, i.group(),i.span()[0], i.span()[1], record.seq.tostring()[i.span()[0]-8:i.span()[0]], record.seq.tostring(), Read_type])
				outfile.writerow([record.id, i.group(),i.span()[0], i.span()[1], record.seq.tostring()[i.span()[0]-8:i.span()[0]], record.seq.tostring(), Read_type])
				SeqIO.write(SeqRecord_gernerator("IGG_"+record.id, record.seq.tostring()[i.span()[0]-8:], Read_type), cut_berfore_UMI_fasta, "fasta")
			if index == 0 and i.span()[0]-8 < 0:
				#print type(i.group()),i.span()
				outfile1.writerow([record.id,i.group(),i.span()[0],i.span()[1],record.seq.tostring()[:i.span()[0]],record.seq.tostring(), Read_type])
				outfile.writerow([record.id,i.group(),i.span()[0],i.span()[1],record.seq.tostring()[:i.span()[0]],record.seq.tostring(), Read_type])
				SeqIO.write(SeqRecord_gernerator("IGG_"+record.id, record.seq.tostring()[:], Read_type), cut_berfore_UMI_fasta, "fasta")
		for index, i in enumerate(loc3):
			if index == 0 and i.span()[0]-8 >= 0:
				#print type(i.group()),i.span()
				outfile3.writerow([record.id,i.group(),i.span()[0],i.span()[1],record.seq.tostring()[i.span()[0]-8:i.span()[0]],record.seq.tostring(), Read_type])
				outfile.writerow([record.id,i.group(),i.span()[0],i.span()[1],record.seq.tostring()[i.span()[0]-8:i.span()[0]],record.seq.tostring(), Read_type])
				SeqIO.write(SeqRecord_gernerator("IGM_"+record.id, record.seq.tostring()[i.span()[0]-8:], Read_type), cut_berfore_UMI_fasta, "fasta")
			if index == 0 and i.span()[0]-8 < 0:
				#print type(i.group()),i.span()
				outfile3.writerow([record.id,i.group(),i.span()[0],i.span()[1],record.seq.tostring()[:i.span()[0]],record.seq.tostring(), Read_type])
				outfile.writerow([record.id,i.group(),i.span()[0],i.span()[1],record.seq.tostring()[:i.span()[0]],record.seq.tostring(), Read_type])
				SeqIO.write(SeqRecord_gernerator("IGM_"+record.id, record.seq.tostring()[:], Read_type), cut_berfore_UMI_fasta, "fasta")
		for index, i in enumerate(loc2):
			if index == 0  and i.span()[0]-8 >= 0:
				#print type(i.group()),i.span()
				outfile2.writerow([record.id,i.group(),i.span()[0],i.span()[1],record.seq.tostring()[i.span()[0]-8:i.span()[0]],record.seq.tostring(), Read_type])
				outfile.writerow([record.id,i.group(),i.span()[0],i.span()[1],record.seq.tostring()[i.span()[0]-8:i.span()[0]],record.seq.tostring(), Read_type])
				SeqIO.write(SeqRecord_gernerator("IGK_"+record.id, record.seq.tostring()[i.span()[0]-8:], Read_type), cut_berfore_UMI_fasta, "fasta")
			if index == 0 and i.span()[0]-8 < 0:
				#print type(i.group()),i.span()
				outfile2.writerow([record.id,i.group(),i.span()[0],i.span()[1],record.seq.tostring()[:i.span()[0]],record.seq.tostring(), Read_type])
				outfile.writerow([record.id,i.group(),i.span()[0],i.span()[1],record.seq.tostring()[:i.span()[0]],record.seq.tostring(), Read_type])	
				SeqIO.write(SeqRecord_gernerator("IGK_"+record.id, record.seq.tostring()[:], Read_type), cut_berfore_UMI_fasta, "fasta")
		for index, i in enumerate(loc4):
			if index == 0 and i.span()[0]-8 >= 0:
				#print type(i.group()),i.span()
				outfile4.writerow([record.id,i.group(),i.span()[0],i.span()[1],record.seq.tostring()[i.span()[0]-8:i.span()[0]],record.seq.tostring(), Read_type])
				outfile.writerow([record.id,i.group(),i.span()[0],i.span()[1],record.seq.tostring()[i.span()[0]-8:i.span()[0]],record.seq.tostring(), Read_type])
				SeqIO.write(SeqRecord_gernerator("IGL_"+record.id, record.seq.tostring()[i.span()[0]-8:], Read_type), cut_berfore_UMI_fasta, "fasta")
			if index == 0 and i.span()[0]-8 < 0:
				#print type(i.group()),i.span()
				outfile4.writerow([record.id,i.group(),i.span()[0],i.span()[1],record.seq.tostring()[:i.span()[0]],record.seq.tostring(), Read_type])
				outfile.writerow([record.id,i.group(),i.span()[0],i.span()[1],record.seq.tostring()[:i.span()[0]],record.seq.tostring(), Read_type])
				SeqIO.write(SeqRecord_gernerator("IGL_"+record.id, record.seq.tostring()[:], Read_type), cut_berfore_UMI_fasta, "fasta")
	
		
def caculate_UMI(fname):
	infile = csv.reader(open("%s"%fname,'rU'),delimiter = '\t')
	outputfile1 = csv.writer(open("%s-count_UMI.txt"%fname,'w'),delimiter = '\t')
	outputfile2 = csv.writer(open("%s-reads_UMI.txt"%fname,'w'),delimiter = '\t')
	outputfile3 = csv.writer(open("%s-reads_no_UMI.txt"%fname,'w'),delimiter = '\t')
	result_dict = {}
	
	for line in infile:
		if "A" in line[4] or "T" in line[4] or "G" in line[4] or "C" in line[4]:
		#if len(line[4]) == 8:
			result_dict.setdefault(str(line[4]), []).append(str(line[0]))
		else:
			
			outputfile3.writerow(line)
	#print result_dict
	for index,key in enumerate(sorted(result_dict.keys())):
		#print key, result_dict[key]
		outputfile1.writerow([key, len(result_dict[key])])
		new_value = copy.deepcopy(result_dict[key])
		new_value.insert(0,key)
		outputfile2.writerow(new_value)
		
def caculate_UMI_8(fname):
	infile = csv.reader(open("%s"%fname,'rU'),delimiter = '\t')
	outputfile1 = csv.writer(open("%s-count_UMI_8.txt"%fname,'w'),delimiter = '\t')
	outputfile2 = csv.writer(open("%s-reads_UMI_8.txt"%fname,'w'),delimiter = '\t')
	outputfile3 = csv.writer(open("%s-reads_no_8_UMI.txt"%fname,'w'),delimiter = '\t')
	result_dict = {}
	
	for line in infile:
		#if "A" in line[4] or "T" in line[4] or "G" in line[4] or "C" in line[4]:
		if len(line[4]) == 8 and ("A" in line[4] or "T" in line[4] or "G" in line[4] or "C" in line[4]):
			result_dict.setdefault(str(line[4]), []).append(str(line[0]))
		else:
			
			outputfile3.writerow(line)
	#print result_dict
	for index,key in enumerate(sorted(result_dict.keys())):
		#print key, result_dict[key]
		outputfile1.writerow([key, len(result_dict[key])])
		new_value = copy.deepcopy(result_dict[key])
		new_value.insert(0,key)
		outputfile2.writerow(new_value)

def caculate_UMI_from6_8(fname):
	infile = csv.reader(open("%s"%fname,'rU'),delimiter = '\t')
	outputfile1 = csv.writer(open("%s-count_UMI_6_8.txt"%fname,'w'),delimiter = '\t')
	outputfile2 = csv.writer(open("%s-reads_UMI_6_8.txt"%fname,'w'),delimiter = '\t')
	outputfile3 = csv.writer(open("%s-reads_no_6_8_UMI.txt"%fname,'w'),delimiter = '\t')
	result_dict = {}
	
	for line in infile:
		#if "A" in line[4] or "T" in line[4] or "G" in line[4] or "C" in line[4]:
		if len(line[4]) >= 6 and len(line[4]) <= 8 and ("A" in line[4] or "T" in line[4] or "G" in line[4] or "C" in line[4]):
			result_dict.setdefault(str(line[4]), []).append(str(line[0]))
		else:
			
			outputfile3.writerow(line)
	#print result_dict
	for index,key in enumerate(sorted(result_dict.keys())):
		#print key, result_dict[key]
		outputfile1.writerow([key, len(result_dict[key])])
		new_value = copy.deepcopy(result_dict[key])
		new_value.insert(0,key)
		outputfile2.writerow(new_value)
def bsub_IgBLAST_jobs(IgBLAST_job):
	cmd = '%s <  %s'%('bsub', IgBLAST_job)
	head, tail 	= os.path.splitext(IgBLAST_job)
	f_ind 		= head.split("_")[ -1 ]
	
	p=subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)  
	while True:
		buff = p.stdout.readline()
		if ">" in buff and "<" in buff:
	 		igblast_job_id = buff.split(">")[0].split("<")[-1]
		if buff == '' and p.poll() != None:
			break
	return igblast_job_id
def check_IgBLAST_done(igblast_job_ids):
	log_file = "%s/IgBLAST_Done.log"%(prj_tree.logs)
	log_file_handle = open(log_file, "w")
	for igblast_job_id in igblast_job_ids:
		while os.path.exists("%s/output_%s"%(prj_tree.home, igblast_job_id)) == False:
			log_file_handle.write("Waiting for IgBLAST_job_%s...%s\n"%(igblast_job_id,time.ctime()))
			time.sleep(18)
		
		errput, output = "%s/errput_%s"%(prj_tree.home, igblast_job_id), "%s/output_%s"%(prj_tree.home, igblast_job_id)
		IgBLAST_log = True
		while IgBLAST_log:
			output_log = open(output, "rU")
			for line in output_log.readlines():
				if line.replace('\n','') == "Successfully completed.":
					log_file_handle.write("IgBLAST_%s Successfully completed.\n"%igblast_job_id)
					IgBLAST_log = False
					os.system("mv %s %s"%(errput, prj_tree.logs))
					os.system("mv %s %s"%(output, prj_tree.logs))
			output_log.close()
def get_assignment_and_recombanation_info(infile):
	fname = retrieve_name_body(infile)
	count_v,count_d,count_j = 0,0,0 
	c_v, c_d, c_j, result ='', '', '', []
	outfile = open("%s/%s_get_assignment_info.txt"%(prj_tree.igblast_data, fname),"w")
	writer = csv.writer(outfile,delimiter = "\t")
	outfile2 = open("%s/%s_get_recombanation_info.txt"%(prj_tree.igblast_data, fname),"w")
	writer2 = csv.writer(outfile2,delimiter = "\t")
	reader = csv.reader(open(infile,"rU"),delimiter = "\t")
	for line in reader:
		con = str(line)
		con=con.replace('\'','')
		con=con[1:-1]
		hit = re.findall('# Query',con)
		hit_v=re.match(r'^V.+:.+',con)
		hit_d=re.match(r'^D.+:.+',con)
		hit_j=re.match(r'^J.+:.+',con)
		if hit_v:
			a_v=hit_v.group()
			b_v = a_v.split(',')
			if not(c_v == b_v[1]):
				c_v=b_v[1]
				count_v += 1
				writer.writerow(b_v)
		if hit_d:
			a_d=hit_d.group()
			b_d = a_d.split(',')
			if not(c_d == b_d[1]):
				c_d=b_d[1]
				count_d += 1
				writer.writerow(b_d)
		if hit_j:
			a_j=hit_j.group()
			b_j = a_j.split(',')
			if not(c_j == b_j[1]):
				c_j=b_j[1]
				count_j += 1
				writer.writerow(b_j)
		
		if hit:
			con = [':'.join(con.split(":")[1:])]
			con[0].replace(" ", "")
			result.append(con)
		if len(line) >= 7:
			if line[-1] == '+' or line[-1] == '-':
				result[-1].extend(line)
	writer2.writerows(result)
	
if __name__=='__main__':
	print 'Parent process %s'%os.getpid()
	
	#step 0
	#rename
	# create 1st and 2nd subfolders
	prj_folder = os.getcwd()
	#prj_tree = create_folders(prj_folder)
	prj_tree = ProjectFolders(prj_folder)
	
	prj_name = fullpath2last_folder(prj_tree.home)
	"""
	#Step 1	Convert fastq to fasta
	print "Step 1	Convert fastq to fasta"
	infile = glob.glob("%s/*.fastq"%(prj_tree.origin))[0]
	fname, suffix = os.path.splitext(infile)
	the_count = SeqIO.convert(infile, 'fastq', '%s/%s.fasta'%(prj_tree.origin, prj_name), "fasta")
	print "Converted %i records to fasta" % (the_count)
	
	#Step 2: Split to little files
	print "Step 2: Split to little files"
	record_iter = SeqIO.parse(open("%s/%s.fasta"%(prj_tree.origin, prj_name)), "fasta")
	for i, batch in enumerate(batch_iterator(record_iter, 50000)) :
		filename = "%s/%s_%i.fasta" % (prj_tree.split,prj_name, i+1)
		handle = open(filename, "w")
		count = SeqIO.write(batch, handle, "fasta")
		handle.close()
		print "Wrote %i records to %s" % (count, filename)
	files_num = i+1
	"""
	'''
	#Step 5: Mapping, Multiple processing
	print "Begin IgBLAST..."
	prepare_IgBLAST_jobs(prj_name, prj_tree)
	IgBLAST_jobs = glob.glob("%s/IgBLAST_*.sh" %(prj_tree.jobs))
	IgBLAST_pool = Pool()
	#for IgBLAST_job in IgBLAST_jobs:
	#	print "IgBLAST_job: %s has been submited."%IgBLAST_job
	IgBLAST_jobs_ids = IgBLAST_pool.map_async(bsub_IgBLAST_jobs, IgBLAST_jobs).get(120)
	print "IgBLAST_jobs: %s IgBLAST_job has been submited."%len(IgBLAST_jobs_ids)
	check_IgBLAST_done(IgBLAST_jobs_ids)
	print "Waiting for all subprocesses done..."
	IgBLAST_pool.close()
	IgBLAST_pool.join()
	print 'All subprocesses done.'
	'''
	
	"""
	pool_size = multiprocessing.cpu_count()
	os.system("rm %s/IgBLAST_result_*get_assignment_info.txt"%(prj_tree.igblast_data))
	os.system("rm %s/IgBLAST_result_*get_recombanation_info.txt"%(prj_tree.igblast_data))
	igblast_result_files = glob.glob("%s/IgBLAST_result_*.txt"%(prj_tree.igblast_data))
	pool = Pool(processes = pool_size-1)
	for igblast_result_file in igblast_result_files:
		
		#get_assignment_and_recombanation_info(igblast_result_file)
		pool.apply_async(get_assignment_and_recombanation_info, args=(igblast_result_file,))
	print "Waiting for all subprocesses done..."
	pool.close()
	pool.join()
	print 'All subprocesses done.'
	os.system("cat %s/IgBLAST_result_*_get_assignment_info.txt > %s/%s_get_assignment_info.txt"%(prj_tree.igblast_data, prj_tree.igblast_data, prj_name))
	os.system("cat %s/IgBLAST_result_*_get_recombanation_info.txt > %s/%s_get_recombanation_info.txt"%(prj_tree.igblast_data, prj_tree.igblast_data, prj_name))
	
	"""
	#files_num = 128
	#for f_idex in range(1, files_num+1):
	#	igblast_result_file = "%s/IgBLAST_result_%s.txt"%(prj_tree.igblast_data,f_idex)
	
	#"""
	#Step 3: detect UMI, split file across UMI 
	pool_size = multiprocessing.cpu_count()
	p1 = Pool(processes = pool_size-1)
	manager = Manager()
	#semaphore = Semaphore(16)
	#print help(manager.dict())
	#sys.exit(0)
	recombanation_info_file = "%s/%s_get_recombanation_info.txt"%(prj_tree.igblast_data, prj_name)
	recombanation_info_dict = load_recombanation_info_dict(recombanation_info_file)
	shared_dict = manager.dict(recombanation_info_dict)
	infiles = glob.glob("%s/%s*.fasta"%(prj_tree.split, prj_name))
	for handle in infiles:
		fname = retrieve_name_body(handle)
		print "Searching...%s"%fname
		#get_primer(handle)
		#get_primer_and_cut_before_UMI(handle, recombanation_info_dict)
		p1.apply_async(get_primer_and_cut_before_UMI, args=(handle, shared_dict))

	print "Waiting for all subprocesses done..."
	p1.close()
	p1.join()
	os.system("cat %s/%s*_cut_before_UMI.fasta >  %s/%s_cut_before_UMI.fasta"%(prj_tree.tmp, prj_name, prj_tree.origin, prj_name))
	print 'Detect primer: All subprocesses done.'
	#"""
	
	'''
	#Step 4
	infiles = glob.glob("%s/%s*-get-primer"%(prj_tree.data,prj_name))
	pool_size = multiprocessing.cpu_count()
	p3 = Pool(processes = pool_size-1)
	for handle in infiles:
		fname, suffix = os.path.splitext(handle)
		#caculate_UMI(fname)
		p3.apply_async(caculate_UMI, args=(fname,))
		p3.apply_async(caculate_UMI_8, args=(fname,))
		p3.apply_async(caculate_UMI_from6_8, args=(fname,))
	print "Waiting for all subprocesses done..."
	p3.close()
	p3.join()
	print 'Detect UMI: All subprocesses done.'
	
	'''
	"""
	# Plot UMI_group_frequency
	UMI_lengths = ["","_6_8","_8"]
	for UMI_length in UMI_lengths:
		UMI_count_files = glob.glob("%s/%s_*-get-primer-count_UMI%s.txt"%(prj_tree.data, prj_name, UMI_length))
		UMI_group_dict, total_reads_number, size_list = {}, 0, []
		for UMI_count_file in UMI_count_files:
			print "PROCESSING %s..."%UMI_count_file
			UMI_count_handle = csv.reader(open(UMI_count_file, "rU"), delimiter="\t")
			for line in UMI_count_handle:
				total_reads_number += int(line[1])
				UMI_group_dict.setdefault(line[0], []).append(int(line[1]))
		print "Length of UMI_group_dict: %s"%len(UMI_group_dict)
		UMI_unqiue_count_file = csv.writer(open("%s/%s_get_primer_count_UMI%s.txt"%(prj_tree.data, prj_name, UMI_length), "w"), delimiter="\t")
		for key, value in UMI_group_dict.items():
			UMI_unqiue_count_file.writerow([key, sum(value)])
			size_list.append(sum(value))
		print "size_list:",size_list
		UMI_size_frequency = number_statistics(size_list)
		print "UMI_size_frequency:",sorted(UMI_size_frequency)
		print "Drawinng ..."
		UMI_size_frequency_file = csv.writer(open("%s/%s_UMI%s_size_frequency.txt"%(prj_tree.data, prj_name, UMI_length), "w"), delimiter="\t")
		x_array, y_array = [], []
		for item in sorted(UMI_size_frequency):
			UMI_size_frequency_file.writerow([item[0],float(item[1])/float(total_reads_number)*100])
			if item[0] <= 2000 and item[1] >=3:
				x_array.append(item[0])
				y_array.append(float(item[1])/float(total_reads_number)*100)
		fig = plt.figure(figsize=(5, 5), dpi=300)
		ax 	= fig.add_subplot(111)
		ax.plot(x_array, y_array, 'o-', markersize=1)
		plt.setp(ax.get_xticklabels(), fontsize=4)	#rotation='vertical', 
		plt.setp(ax.get_yticklabels(), fontsize=4)
		ax = plt.gca()
		for I in ax.get_xticklines() + ax.get_yticklines():
			I.set_markersize(2)
		plt.savefig("%s/%s_UMI%s_group_frequency.png"%(prj_tree.figure, prj_name, UMI_length), dpi=300)
	"""
	
	
	
