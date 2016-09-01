#!/usr/bin/env python
# encoding: utf-8
"""
common_info.py

Created by Zhenhai Zhang on 2011-04-06.
Copyright (c) 2011 Columbia University and Vaccine Research Center, National Institutes of Health, USA. All rights reserved.
"""
UMI_length = int(8)

IGG_primer = 'AAGACCGATGGGCCCTTGGTGGA'
IGK_primer = 'AAGACAGATGGTGCAGCCACAGTTC'
IGM_primer = 'AAGGGTTGGGGCGGATGCACTCC'
IGL_primer = 'CCTTGTTGGCTTG(A|G)AGCTCCTCAGAGGAGG'


primer_len_dict = {
	"IGG" : int(len(IGG_primer)),
	"IGK" : int(len(IGK_primer)),
	"IGM" : int(len(IGM_primer)),
	"IGL" : int(len(IGL_primer)-4),
	"No_type" : int(0),
	}


ORG_FOLDER = "origin"
FILTERED_FOLDER="filtered"
MAPPING_FOLDER="mapping" 
ANALYSIS_FOLDER="analysis"
LOG_FOLDER="logs"
CLUSTAL_FOLDER="clustal"
PHYLO_FOLDER="phylo"
DOC_FOLDER="docs"
TMP_FOLDER="tmp"

SPLIT_FOLDER="split" 
GERM_FOLDER="germ"
NAT_FOLDER="native"	
SELF_FOLDER="reads"
PBS_FOLDER="pbs"
JOB_FOLDER="jobs"

ANALYSIS_DATA_FOLDER="data"
ANALYSIS_FIGURE_FOLDER="figure"

CLUSTAL_FASTA_FOLDER="clustal_fasta"
CLUSTAL_PBS_FOLDER="clustal_pbs"
CLUSTAL_JOB_FOLDER="clustal_job"
CLUSTAL_DATA_FOLDER ="clustal_data"

IGBLAST_DATA_FOLDER = "Igblast_data"
IGBLAST_DATABASE_FOLDER = "Igblast_database"
IGBLAST_TMP_FOLDER = "Igblast_tmp"


nucleotides = ["A", "C", "G", "T"]

#
# ===START=== amino acid
#
dict_codon2aa = {

	"ATT" : ("I", "Isoleucine"),
	"ATA" : ("I", "Isoleucine"),
	"ATC" : ("I", "Isoleucine"),
	
	"CTT" : ("L", "Leucine"),
	"CTC" : ("L", "Leucine"),
	"CTA" : ("L", "Leucine"),
	"CTG" : ("L", "Leucine"),
	"TTA" : ("L", "Leucine"),
	"TTG" : ("L", "Leucine"),
	
	"GTT" : ("V", "Valine"),
	"GTC" : ("V", "Valine"),
	"GTA" : ("V", "Valine"),
	"GTG" : ("V", "Valine"),
	
	"TTT" : ("F", "Phenylalanine"),
	"TTC" : ("F", "Phenylalanine"),
	
	"ATG" : ("M", "Methionine"),
	
	"TGT" : ("C", "Cysteine"),
	"TGC" : ("C", "Cysteine"),
	
	"GCT" : ("A", "Alanine"),
	"GCC" : ("A", "Alanine"),
	"GCA" : ("A", "Alanine"),
	"GCG" : ("A", "Alanine"),
	
	"GGT" : ("G", "Glycine"),
	"GGC" : ("G", "Glycine"),
	"GGA" : ("G", "Glycine"),
	"GGG" : ("G", "Glycine"),
	
	"CCT" : ("P", "Proline"),
	"CCC" : ("P", "Proline"),
	"CCA" : ("P", "Proline"),
	"CCG" : ("P", "Proline"),
	
	"ACT" : ("T", "Threonine"),
	"ACC" : ("T", "Threonine"),
	"ACA" : ("T", "Threonine"),
	"ACG" : ("T", "Threonine"),
	
	"TCT" : ("S", "Serine"),
	"TCC" : ("S", "Serine"),
	"TCA" : ("S", "Serine"),
	"TCG" : ("S", "Serine"),
	"AGT" : ("S", "Serine"),
	"AGC" : ("S", "Serine"),
	
	"TAT" : ("Y", "Tyrosine"),
	"TAC" : ("Y", "Tyrosine"),
	
	"TGG" : ("W", "Tryptophan"),
	
	"CAA" : ("Q", "Glutamine"),
	"CAG" : ("Q", "Glutamine"),
	
	"AAT" : ("N", "Asparagine"), 
	"AAC" : ("N", "Asparagine"),
	
	"CAT" : ("H", "Histidine"), 
	"CAC" : ("H", "Histidine"),
	
	"GAA" : ("E", "Glutamic acid"),
	"GAG" : ("E", "Glutamic acid"),
	
	"GAT" : ("D", "Aspartic acid"),
	"GAC" : ("D", "Aspartic acid"),
	
	"AAA" : ("K", "Lysine"),
	"AAG" : ("K", "Lysine"),
	
	"CGT" : ("R", "Arginine"),
	"CGC" : ("R", "Arginine"),
	"CGA" : ("R", "Arginine"),
	"CGG" : ("R", "Arginine"),
	"AGA" : ("R", "Arginine"),
	"AGG" : ("R", "Arginine"),
	
	"TAA" : ("STOP", "Stop codons"),
	"TAG" : ("STOP", "Stop codons"),
	"TGA" : ("STOP", "Stop codons")

	}
	
start_codons = {"ATG"}
stop_codons = {"TAA", "TAG", "TGA"}

dict_aa2codon = {

	"I" : ("ATT", "ATC", "ATA"),
	"L" : ("CTT", "CTC", "CTA", "CTG", "TTA", "TTG"),
	"V" : ("GTT", "GTC", "GTA", "GTG"),
	"F" : ("TTT", "TTC"),
	"M" : ("ATG"),
	"C" : ("TGT", "TGC"),
	"A" : ("GCT", "GCC", "GCA", "GCG"),
	"G" : ("GGT", "GGC", "GGA", "GGG"),
	"P" : ("CCT", "CCC", "CCA", "CCG"),
	"T" : ("ACT", "ACC", "ACG", "ACA"),
	"S" : ("TCT", "TCC", "TCG", "TCA", "AGT", "AGC"),
	"Y" : ("TAT", "TAC"),
	"W" : ("TGG"),
	"Q" : ("CAA", "CAG"),
	"N" : ("AAT", "AAC"),
	"H" : ("CAT", "CAC"),
	"E" : ("GAA", "GAG"),
	"D" : ("GAT", "GAC"),
	"K" : ("AAA", "AAG"),
	"R" : ("CGT", "CGC", "CGA", "CGG", "AGA", "AGG"),
	"STOP" : ("TAA", "TAG", "TGA")

	}
dict_aa2name = {

	"I" : "Isoleucine",
	"L" : "Leucine",
	"V" : "Valine",
	"F" : "Phenylalanine",
	"M" : "Methionine",
	"C" : "Cysteine",
	"A" : "Alanine",
	"G" : "Glycine",
	"P" : "Proline",
	"T" : "Threonine",
	"S" : "Serine",
	"Y" : "Tyrosine",
	"W" : "Tryptophan",
	"Q" : "Glutamine",
	"N" : "Asparagine",
	"H" : "Histidine",
	"E" : "Glutamic acid",
	"D" : "Aspartic acid",
	"K" : "Lysine",
	"R" : "Arginine",
	"STOP" : "Stop condons"
}

#
# ===END=== amino acid
#


#
# ===START=== amino acid
#
IUB_CODE = {
	"A" : ["A"], 
	"C" : ["C"], 
	"G" : ["G"], 
	"T" : ["T"], 
	"R" : ["A", "G"], 
	"Y" : ["C", "T"], 
	"K" : ["G", "T"], 
	"M" : ["A", "C"], 
	"S" : ["G", "C"], 
	"W" : ["A", "T"], 
	"B" : ["C", "G", "T"], 
	"D" : ["A", "G", "T"], 
	"H" : ["A", "C", "T"], 
	"V" : ["A", "C", "G"], 
	"N" : ["A", "C", "G", "T"]
}
	
IUB_complement = {
	"A" : "T",
	"C" : "G", 
	"G" : "C", 
	"T" : "A",
	"R" : "Y", 
	"Y" : "R",
	"K" : "M", 
	"M" : "K",
	"S" : "W",
	"W" : "S",
	"B" : "V",
	"V" : "B",
	"D" : "H",
	"H" : "D",
	"N" : "N"
}
#
# ===END=== amino acid
#
