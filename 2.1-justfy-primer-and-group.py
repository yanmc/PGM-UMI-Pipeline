#!/usr/bin/env python
# encoding: utf-8
"""
2.1_justfy_primer_and_group.py

Created by Mingchen on 2015-01-19.
Copyright (c) 2015 __MyCompanyName__. All rights reserved.

Usage: misc-get-pattern-loc.py -i infile -r referencefile -o outfile


"""
def justfy_primer_and_group(key,value,origin_fasta_dict, recombanation_info_dict, UMI_length):
	#print 1
	print "Write No.%s UMI, it has %s reads."%(index, len(value))
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

if __name__=='__main__':
	print 'Parent process %s'%os.getpid()
	prj_folder = os.getcwd()
	
	prj_tree = ProjectFolders(prj_folder)
	
	prj_name = fullpath2last_folder(prj_tree.home)
	dict_args = processParas(sys.argv, i="infile")
	infile = getParas(dict_args, "infile")
	#'''
	# Step1: Split same UMI 8 reads to same file
	print "Loading origin fasta dict......"
	
	#origin_fasta_dict = SeqIO.index("%s/%s_cut_before_UMI.fasta"%(prj_tree.origin, prj_name), "fasta")
	origin_fasta_dict = load_fastas("%s/%s_cut_before_UMI.fasta"%(prj_tree.origin, prj_name))
	print "reads number:%s"%(len(origin_fasta_dict))
	
	recombanation_info_file = "%s/%s_get_recombanation_info.txt"%(prj_tree.igblast_data, prj_name)
	recombanation_info_dict = load_recombanation_info_dict(recombanation_info_file)
	f = open(infile, 'rb')
	pickle_tuple = pickle.load(f)
	key, value = pickle_tuple[0], pickle_tuple[1]
	f.close()
	
	#UMI_lengths = ["","_6_8","_8"]
	UMI_lengths = ["_8"]
	for UMI_length in UMI_lengths:
			justfy_primer_and_group(key,value,origin_fasta_dict, recombanation_info_dict, UMI_length)
	print "Done!"
			#group_pool.apply_async(justfy_primer_and_group, args= (key, value, shared_origin_fasta_dict, shared_recombanation_info_dict, UMI_length))