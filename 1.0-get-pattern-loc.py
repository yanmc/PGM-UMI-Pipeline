#!/usr/bin/env python
# encoding: utf-8
"""
misc-get-pattern-loc.py

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
def get_primer(handle):
	file_handle = open("%s"%handle, 'rU')
	fname = retrieve_name_body(handle)
	outfile  = csv.writer(open("%s/%s-get-primer"%(prj_tree.data,fname),'w'),delimiter = '\t')
	outfile1 = csv.writer(open("%s/%s-get-IgG-primer"%(prj_tree.tmp,fname),'w'),delimiter = '\t')
	outfile2 = csv.writer(open("%s/%s-get-IgK-primer"%(prj_tree.tmp,fname),'w'),delimiter = '\t')
	outfile3 = csv.writer(open("%s/%s-get-IgM-primer"%(prj_tree.tmp,fname),'w'),delimiter = '\t')
	outfile4 = csv.writer(open("%s/%s-get-IgL-primer"%(prj_tree.tmp,fname),'w'),delimiter = '\t')
	#outfile1.writerow(['id','primer','start','end','UMI',"Seq"])
	#outfile2.writerow(['id','primer','start','end','UMI',"Seq"])
	#outfile3.writerow(['id','primer','start','end','UMI',"Seq"])
	#outfile4.writerow(['id','primer','start','end','UMI',"Seq"])
	#print "Searching...%s"%fname
	the_handle = SeqIO.parse(file_handle, 'fasta')
	for record in the_handle:
		#print "Searching...%s"%record.id
		#print record.seq.tostring()
		#loc = re.finditer('AAGACCGATGGGCCCTTGGTGGA',record.seq.reverse_complement().tostring(),re.I)
		loc1 = re.finditer('AAGACCGATGGGCCCTTGGTGGA',record.seq.tostring(),re.I)
		loc2 = re.finditer('AAGACAGATGGTGCAGCCACAGTTC',record.seq.tostring(),re.I)
		loc3 = re.finditer('AAGGGTTGGGGCGGATGCACTCC',record.seq.tostring(),re.I)
		loc4 = re.finditer('CCTTGTTGGCTTG(A|G)AGCTCCTCAGAGGAGG',record.seq.tostring(),re.I)
		#loc = re.finditer('cacagtg',record.seq.tostring())
		#loc = re.finditer('ACAAAAACC',record.seq.tostring(),re.I)
		for index, i in enumerate(loc1):
			if index == 0 and i.span()[0]-8 >= 0:
				#print type(i.group()),i.span()
				outfile1.writerow([record.id,i.group(),i.span()[0],i.span()[1],record.seq.tostring()[i.span()[0]-8:i.span()[0]],record.seq.tostring()])
				outfile.writerow([record.id,i.group(),i.span()[0],i.span()[1],record.seq.tostring()[i.span()[0]-8:i.span()[0]],record.seq.tostring()])
			if index == 0 and i.span()[0]-8 < 0:
				#print type(i.group()),i.span()
				outfile1.writerow([record.id,i.group(),i.span()[0],i.span()[1],record.seq.tostring()[:i.span()[0]],record.seq.tostring()])
				outfile.writerow([record.id,i.group(),i.span()[0],i.span()[1],record.seq.tostring()[:i.span()[0]],record.seq.tostring()])
		for index, i in enumerate(loc2):
			if index == 0  and i.span()[0]-8 >= 0:
				#print type(i.group()),i.span()
				outfile2.writerow([record.id,i.group(),i.span()[0],i.span()[1],record.seq.tostring()[i.span()[0]-8:i.span()[0]],record.seq.tostring()])
				outfile.writerow([record.id,i.group(),i.span()[0],i.span()[1],record.seq.tostring()[i.span()[0]-8:i.span()[0]],record.seq.tostring()])
			if index == 0 and i.span()[0]-8 < 0:
				#print type(i.group()),i.span()
				outfile2.writerow([record.id,i.group(),i.span()[0],i.span()[1],record.seq.tostring()[:i.span()[0]],record.seq.tostring()])
				outfile.writerow([record.id,i.group(),i.span()[0],i.span()[1],record.seq.tostring()[:i.span()[0]],record.seq.tostring()])	
		for index, i in enumerate(loc3):
			if index == 0 and i.span()[0]-8 >= 0:
				#print type(i.group()),i.span()
				outfile3.writerow([record.id,i.group(),i.span()[0],i.span()[1],record.seq.tostring()[i.span()[0]-8:i.span()[0]],record.seq.tostring()])
				outfile.writerow([record.id,i.group(),i.span()[0],i.span()[1],record.seq.tostring()[i.span()[0]-8:i.span()[0]],record.seq.tostring()])
			if index == 0 and i.span()[0]-8 < 0:
				#print type(i.group()),i.span()
				outfile3.writerow([record.id,i.group(),i.span()[0],i.span()[1],record.seq.tostring()[:i.span()[0]],record.seq.tostring()])
				outfile.writerow([record.id,i.group(),i.span()[0],i.span()[1],record.seq.tostring()[:i.span()[0]],record.seq.tostring()])
		for index, i in enumerate(loc4):
			if index == 0 and i.span()[0]-8 >= 0:
				#print type(i.group()),i.span()
				outfile4.writerow([record.id,i.group(),i.span()[0],i.span()[1],record.seq.tostring()[i.span()[0]-8:i.span()[0]],record.seq.tostring()])
				outfile.writerow([record.id,i.group(),i.span()[0],i.span()[1],record.seq.tostring()[i.span()[0]-8:i.span()[0]],record.seq.tostring()])
			if index == 0 and i.span()[0]-8 < 0:
				#print type(i.group()),i.span()
				outfile4.writerow([record.id,i.group(),i.span()[0],i.span()[1],record.seq.tostring()[:i.span()[0]],record.seq.tostring()])
				outfile.writerow([record.id,i.group(),i.span()[0],i.span()[1],record.seq.tostring()[:i.span()[0]],record.seq.tostring()])

def caculate_UMI(fname):
	infile = csv.reader(open("%s.txt"%fname,'rU'),delimiter = '\t')
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
	infile = csv.reader(open("%s.txt"%fname,'rU'),delimiter = '\t')
	outputfile1 = csv.writer(open("%s-count_UMI_8.txt"%fname,'w'),delimiter = '\t')
	outputfile2 = csv.writer(open("%s-reads_UMI_8.txt"%fname,'w'),delimiter = '\t')
	outputfile3 = csv.writer(open("%s-reads_no_8_UMI.txt"%fname,'w'),delimiter = '\t')
	result_dict = {}
	
	for line in infile:
		#if "A" in line[4] or "T" in line[4] or "G" in line[4] or "C" in line[4]:
		if len(line[4]) == 8:
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
	infile = csv.reader(open("%s.txt"%fname,'rU'),delimiter = '\t')
	outputfile1 = csv.writer(open("%s-count_UMI_6_8.txt"%fname,'w'),delimiter = '\t')
	outputfile2 = csv.writer(open("%s-reads_UMI_6_8.txt"%fname,'w'),delimiter = '\t')
	outputfile3 = csv.writer(open("%s-reads_no_6_8_UMI.txt"%fname,'w'),delimiter = '\t')
	result_dict = {}
	
	for line in infile:
		#if "A" in line[4] or "T" in line[4] or "G" in line[4] or "C" in line[4]:
		if len(line[4]) >= 6 and len(line[4]) <= 8 :
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
if __name__=='__main__':
	print 'Parent process %s'%os.getpid()
	
	#step 0
	#rename
	# create 1st and 2nd subfolders
	prj_folder = os.getcwd()
	prj_tree = create_folders(prj_folder)
	#folder_tree = ProjectFolders(prj_folder)
	
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
	
	#Step 3: detect UMI, split file across UMI 
	
	p1 = Pool()
	infiles = glob.glob("%s/%s*.fasta"%(prj_tree.split, prj_name))
	for handle in infiles:
		fname = retrieve_name_body(handle)
		print "Searching...%s"%fname
		p1.apply_async(get_primer, args=(handle,))

	print "Waiting for all subprocesses done..."
	p1.close()
	p1.join()
	print 'Detect primer: All subprocesses done.'
	
	#Step 4
	infiles = glob.glob("%s/%s-get-primer"%(prj_tree.data,fname))
	p3 = Pool()
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
	"""
	
	
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
	#Step 4 analysis IgBLAST result in prj_tree.igblast_data
	IgBLAST_result_files = glob.glob("%s/IgBLAST_result_*"%prj_tree.igblast_data)
	for IgBLAST_result_file in IgBLAST_result_files:
		handle = csv.reader(open(IgBLAST_result_file, "rU"), delimiter='\t')
		for index, line in enumerate(handle):
			print line
			if line != [] and "# Query:" in line[0]:
				line_id_marker = index
				V_alignment, D_alignment, J_alignment = '', '', ''

			if line != [] and line[0] == "V":
				V_alignment = MyAlignment(line)
			if line != [] and line[0] == "D":
				D_alignment = MyAlignment(line)
			if line != [] and line[0] == "J":
				J_alignment = MyAlignment(line)
	
	
	#Step 2
	p2 = Pool()
	for chain_type in ["IgG","IgK","IgL","IgM"]:
		#combine_file(chain_type)
		p2.apply_async(combine_file, args=(chain_type,))
	
	print "Waiting for all subprocesses done..."
	p2.close()
	p2.join()
	
	#Step 3
	infiles = glob.glob("*-primer")
	p3 = Pool()
	for handle in infiles:
		fname, suffix = os.path.splitext(handle)
		#caculate_UMI(fname)
		p3.apply_async(caculate_UMI, args=(fname,))
		p3.apply_async(caculate_UMI_8, args=(fname,))
		p3.apply_async(caculate_UMI_from6_8, args=(fname,))
	print "Waiting for all subprocesses done..."
	p3.close()
	p3.join()
	print 'All subprocesses done.'
	'''
	
