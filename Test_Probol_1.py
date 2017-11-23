#!/usr/bin/env python
# -*- coding:UTF-8 -*-

__author__ = "Zongyun Qiao"
__copyright__ = "Copyright 2017, A Biotech"
__credits__ = [
    "Zongyun Qiao"]  # remember to add yourself
__license__ = "GPL"
__version__ = "0.1-dev, 20171123"
__maintainer__ = "Zongyun Qiao"
__email__ = "gulile@yeah.net"


import os
import sys
import subprocess
import shlex
from Bio.Blast import NCBIXML

def retrieve_Mis_pos(HspM, label):
	p_index = -1

	posList = []
	while True:
		p_index = HspM.find(label, p_index + 1)
		if p_index == -1:
			break
		posList.append(p_index)
	if posList == []: 
		return ["-"]

	return posList

def generate_fa_or_dimerSeq(tableFile):
	res_handle = open(tableFile, "r")
	
	x_Table = ".".join(tableFile.split(".")[:-1])
	
	outFiles = (x_Table + ".fulllen.fa", x_Table + ".16bp.fa", x_Table + "_homodimer_In.txt" )
	fullLenFa = open(outFiles[0], "w")
	L16_Fa    = open(outFiles[1], "w")
	HomoDimer_In = open(outFiles[2] ,"w")
	
	for line in res_handle:
		lin_spl_F = line.rstrip().split("\t") 
		fullLenFa.write(">{0}_F\n{1}\n".format(lin_spl_F[0], lin_spl_F[4]))
		L16_Fa.write(">{0}\n{1}\n".format(lin_spl_F[0], lin_spl_F[8]))
		HomoDimer_In.write("{0}\t{1}\n".format(lin_spl_F[0], lin_spl_F[4]))
		
	fullLenFa.close()
	L16_Fa.close()
	HomoDimer_In.close()
	
	return outFiles
	
class FileParser:
	def __init__(self, filename):
		self.f = filename

	def writeHeader(self, input_fmt, outName="test_out.txt" ):
		outF = open(outName, "w")
		if input_fmt == "xml":
			outF.write("""QuerySeq_ID\tquery_length\thits_alignmentID\thsp_matchlabel\tquery_start\
                \tmismatch_Ofalignment\tquery_end\tMismatch_From3end\talignment_length\
                \tright_mismatchNumber\tChromosomeID\tsubject_start\tsubject_end\tquery_seq\tsubject_seq\ 
                \tHSP_bitScore\thsp_score\thsp_eValue\n""")

		outF.close()

	def parse_xml(self, outName):
		out_handle = open(outName, "a")
		
		xmlIn = open(self.f, "r")
		
		blast_records = NCBIXML.parse(xmlIn)
		for record in blast_records:

			for i , alignment in enumerate(record.alignments):
				try:
			## out_align.write(alignment.title + "--------------\n")
					for hsp in alignment.hsps:
						out_handle.write(record.query + "\t"+ str(record.query_length)+"\t" + alignment.title + "\t" )

						out_handle.write(hsp.match + "\t" )
				## out_align.write(str(alignment.length))
				
						if hsp.query_end < hsp.query_start:
							HSPQ_start = hsp.query_end
							HSPQ_end   = hsp.query_start
                      
						else:
							HSPQ_start = hsp.query_start
							HSPQ_end   = hsp.query_end
				
						if hsp.match.find(" ") == -1: ## MisPP is postions of mismatch in the hit query sequence
							misPL = "-"
							right_PL = "="
						else:
							if hsp.query_end < hsp.query_start:

								 revMatch = hsp.match[::-1]
								 misPs = retrieve_Mis_pos(revMatch, " ")                       
							else:

								 misPs = retrieve_Mis_pos(hsp.match, " ")
						 
							misPnum = [i+int(HSPQ_start) for i in misPs]
							misPL = ",".join([str(i+int(HSPQ_start)) for i in misPs])
							right_PL = ",".join([str(record.query_length +1 - i) for i in misPnum])
				
						if alignment.title.find("chromosome") != -1:
							z = alignment.title.find("chromosome")
							chrName = "chr" + alignment.title[z+11:].split()[0].strip(",")
						else:
							chrName = alignment.title
				
						if chrName.find("chr") != -1:
                    
							cIndex = chrName.find("chr")
							chrName = chrName[cIndex:]
                    
							out_handle.write(str(HSPQ_start) +"\t" + misPL + "\t" + str(HSPQ_end)+"\t" + right_PL + "\t" + str(hsp.align_length) + "\t" +str(record.query_length - HSPQ_end ) + "\t" + chrName +"\t" + str(hsp.sbjct_start) + "\t" + str(hsp.sbjct_end) + "\t"  + hsp.query + "\t" + hsp.sbjct + "\t")
				
							out_handle.write("{0}\t{1}\t{2}".format(hsp.bits, hsp.score, hsp.expect))
						out_handle.write("\n")
				except:
					out_handle.write("No Hits\tNone\n")
		out_handle.close()
		xmlIn.close()

	def SaveData(self, record, outName):
		with open(outName, "w") as outF:
			outF.write(record + "\n")

	@staticmethod
	def GetBlastScore(parse_results):
		
		pass
		#return dict
    
	@staticmethod
	def GetDeltaG(dimer_results):
		pass
		#return dict

class ProbeSolution:

	
	def __init__(self, originalSeq, hg19_db = "/data/qiaozongyun/Genome_reference/HG19_BlastDB" ):
		self.seq = originalSeq.upper()
		self.db = hg19_db
		
	def subSeq(self, seqLen, lastNum):
		fullLen = len(self.seq)
		for i in range(0, fullLen - seqLen):
			spliceSeq = self.seq[i:i+seqLen]
			LastBases = spliceSeq[-lastNum:]
			yield (i, i+seqLen, spliceSeq , LastBases)

	def do_Blastn(self, blastProg, seqa_file):
		b_command = shlex.split("{0} -task blastn -query {1} -db {2} -out {1}_blastn.xml -outfmt 5 -max_target_seqs 100".format(blastProg, seqa_file, self.db, seqa_file))
		subprocess.call(b_command)

	def do_Blastall(self, blastProg, seqb_file):
		b_command = shlex.split("{0} -p blastn -i {1} -d {2} -o {1}_blastall.xml -m 7 -b 100".format(blastProg, seqb_file, self.db, seqb_file))
		subprocess.call(b_command)
    
    
	@staticmethod
	def RevComp(seq):
		dBases = {"A":"T", "T":"A", "C":"G", "G":"C"}
		oSeq = ""
		for base in seq:
			oSeq += dBases[base]
		return oSeq[::-1]

	@staticmethod
	def GCCalc(seq):
		GC_Content = seq.count("C") + seq.count("G")
		return GC_Content * 100.0 / len(seq)
    
	@staticmethod
	def TemMelt(seq, minT, maxT):
		sGC = (seq.count("C") + seq.count("G"))*1.0 / len(seq)
		con005_T = 59.9 + 41.0 * sGC - (675.0/len(seq))
		
		con1_T   = 81.5 + 41.0 * sGC - (675.0/len(seq))
		
		if con005_T >= minT * 1.0 and con005_T <= maxT * 1.0:
			return True,  con005_T, con1_T
		else:
			return False, con005_T, con1_T
		 

class Table_Filter:
	def __init__(self , tabFile):
		self.tmp = tabFile
		
	def filter_blast(self, blastScore_dict, bS = 49.5 ):
		pass
		#return list
    
	def filter_selfDimer(self, DeltaG_dict, dG = -5.6):
		pass
		# return list
	def combine_Blastn_dG(self, blastList, dimerList, outFile):
		pass
		# return None

if __name__ == '__main__':

    
	#SeqManipu = ProbeSolution("CCTGGCAGCCAGGAACGTACTGGTGAAAACACCGCAGCATGTCAAGATCACAGATTTTGGGCT")
	SeqManipu = ProbeSolution("""CATGAACTACTTGGAGGACCGTCGCTTGGTGCACCGCGACCTGGCAGCCAGGA\
ACGTACTGGTGAAAACACCGCAGCATGTCAAGATCACAGATTTTGGGC""")   #hg19_dna range=chr7:55259414-55259514 5'pad=0 strand=+
    
	maxLen, minLen = 40, 25
	
	PLabel = "EGFR_PT"
	
	out_test = open("TestEGFR_pt.xls", "w")
	
	for seq_len in range(minLen, maxLen + 1):
		for istart, iend, s, lastS in SeqManipu.subSeq(seq_len, 16):
			InTmRange, Tm005M, Tm1M = ProbeSolution.TemMelt(s, 53, 60)
			if InTmRange:
				#print(ProbeSolution.RevComp(s), Tm005M, Tm1M, ProbeSolution.GCCalc(s))
				out_test.write("{0}_{1}_{2}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\n".format(PLabel, istart, iend, seq_len, s, Tm005M, Tm1M, ProbeSolution.GCCalc(s), lastS, ProbeSolution.GCCalc(lastS)))
	
	out_test.close()
	#generate_fa_or_dimerSeq(tm_GC_table_results)
	generate_fa_or_dimerSeq("TestEGFR_pt.xls")
	
	# do blastn
	SeqManipu.do_Blastall("/thinker/dstore/r3data/qiaozy2/blast-2.2.26/bin/blastall", "TestEGFR_pt.16bp.fa")
	# /home/qiaozy2/Software/ncbi-blast-2.6.0+/bin/blastn -- blastn
	# calculate selfdimer DeltaG   


	#parse blast results
	#ParserA = FileParser("EGFR_blastOut.xml")
    
	#ParserA.writeHeader("xml",  "Class_parser_testOut.xls" )
    
	#ParserA.parse_xml("Class_parser_testOut.xls")

	# get blast score     # how can we define blast score ?????????????????????????????
	# get deltaG
	
	# filter with blast
	
	
	# filter with dG

	#do combine_Blastn_dG evaluation
	
