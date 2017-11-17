#!/usr/bin/env python
# -*- coding:UTF-8 -*-

__author__ = "Zongyun Qiao"
__copyright__ = "Copyright 2017, A Biotech"
__credits__ = [
    "Zongyun Qiao"]  # remember to add yourself
__license__ = "GPL"
__version__ = "0.1-dev, 20171116"
__maintainer__ = "Zongyun Qiao"
__email__ = "gulile@yeah.net"


import os
import sys
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

class ProbeSolution:

	
	def __init__(self, originalSeq ):
		self.seq = originalSeq.upper()
		
	def subSeq(self, seqLen, lastNum):
		fullLen = len(self.seq)
		for i in range(0, fullLen - seqLen):
			spliceSeq = self.seq[i:i+seqLen]
			LastBases = spliceSeq[-lastNum:]
			yield (i, i+seqLen, spliceSeq , LastBases)

	def trim3end_ToBlastn(self):
		pass

	def trim5end_ToBlastn(self):
		pass
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
		 
'''
class Table_Filter:
	def __init__(self ):
		
'''

if __name__ == '__main__':
    #ParserA = FileParser("EGFR_blastOut.xml")
    
    #ParserA.writeHeader("xml",  "Class_parser_testOut.xls" )
    
    #ParserA.parse_xml("Class_parser_testOut.xls")
    
	SeqManipB = ProbeSolution("CCTGGCAGCCAGGAACGTACTGGTGAAAACACCGCAGCATGTCAAGATCACAGATTTTGGGCT")
    
	maxLen, minLen = 40, 25
	
	for seq_len in range(minLen, maxLen + 1):
		for istart, iend, s, lastS in SeqManipB.subSeq(seq_len, 16):
			InTmRange, Tm005M, Tm1M = ProbeSolution.TemMelt(s, 53, 60)
			if InTmRange:
				#print(ProbeSolution.RevComp(s), Tm005M, Tm1M, ProbeSolution.GCCalc(s))
				print(istart, iend, seq_len, s, Tm005M, Tm1M, ProbeSolution.GCCalc(s))
