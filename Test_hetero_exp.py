#!/usr/bin/env python
# -*- coding:UTF-8 -*-

__author__ = "Zongyun Qiao"
__copyright__ = "Copyright 2017, A Biotech"
__credits__ = [
    "Zongyun Qiao"]  # remember to add yourself
__license__ = "GPL"
__version__ = "0.1-dev, 20171206"
__maintainer__ = "Zongyun Qiao"
__email__ = "gulile@yeah.net"


import os
import sys
import subprocess
import shlex
import math
from collections import defaultdict
from collections import Counter
from Bio.Blast import NCBIXML


Q_LA = 0.5967479540423579         # lambda for Expect value calculation
hg19_Len = 3095693983
R_K  = 0.6293523487909248         # K for Expect value calculation

def get_probeSeq(infile):
	d_seqs = {}
	f_handle = open(infile, "r")
	for line in f_handle:
		if line.startswith("gene_name"):
			continue
		else:
			sf = line.rstrip().split("\t")
			seqName = "testP_" + sf[1]
			if sf[2] != "ND":
				d_seqs[seqName] = sf[2]
	f_handle.close()
	return d_seqs

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
	
def calc_PscoreNP(sL, e5, mismatch, un3, startS = 3.5, S5 = 0.5, reward = 1):
	
	# the function derived from the equation: FinalS = reward * (PT of 3' score of each base, startS is first PT of 3' end) - 3 * (mismatch base number) + (5' match base longer than 16bp) * S5
	
	PTs = {}
	sscale = startS ** (1.0/16)

	for i in range(1, 17):
		PTs[i] = startS 
		startS = startS / sscale

	for j in range(17, 50):
		PTs[j] = S5

	s = 0
	if mismatch.find("=") != -1:
		p = []
	else:        
		misPos = mismatch.split(",")
		p = [int(i) for i in misPos]

	end5 = sL - e5 + 1

	for u in range(un3 +1, end5+1):
		if u not in p:
			s += PTs[u] * reward
		else:
			#s -= PTs[u] * 3
			s -= 3

	return s


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
	def GetBlastScore(parse_source, blastS = 16.0):
		parse_results = parse_source + ".P1score.xls"

		file_handle = open(parse_results , "w")
    
		scoreD = defaultdict(list)

		oF = open(parse_source, "r")

		for line in oF:
			if not line.startswith("QuerySeq_"):
				xline   = line.rstrip().split("\t")
            
				seqId = xline[0]
            
				seqLen  = int(xline[1])
				end5    = int(xline[4])
            
				un3end = int(xline[9])
            
				s = calc_PscoreNP(seqLen, end5, xline[7], un3end)
				qe = R_K * seqLen * hg19_Len * math.pow(10, -Q_LA * s)
            
				if s > blastS:
            
					file_handle.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\n".format(xline[0], xline[1], end5, xline[6], xline[7], xline[9],xline[15], xline[16], xline[17], xline[3], s, qe))
            
					scoreD[seqId].append(s)

		file_handle.close()
		oF.close()

		return parse_results, scoreD
    
	@staticmethod
	def GetDeltaG(dimer_results):

		dd = {}                                         ## dd will store deltaG of each comparison
		openF = open(dimer_results, "r")
		for line in openF:
			if not line.startswith("#") and line.strip() != "":
				z = line.strip().split("\t")
				cID = z[0]
				cDeltaG = z[-1].lstrip("deltaG=")
				if cID not in dd or (cID in dd and float(cDeltaG) < dd[cID]):
					dd[cID] = float(cDeltaG)
		openF.close()
                
		return dd

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
		outxml = seqa_file.strip(".fa") + ".xml"
		b_command = shlex.split("{0} -task blastn -query {1} -db {2} -out {3} -outfmt 5 -evalue 2.0 -num_threads 8 -max_hsps 80".format(blastProg, seqa_file, self.db, outxml))
		subprocess.call(b_command)

	def do_Blastall(self, blastProg, seqb_file):
		outxml = seqb_file.strip(".fa") + "_blastall.xml"
		b_command = shlex.split("{0} -p blastn -i {1} -d {2} -o {3} -a 4 -W 13 -e 2.0 -m 7 -K 1".format(blastProg, seqb_file, self.db, outxml))
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
		
	def filter_blast( self, blastScore, hitNum = 3, ScoreDrop = 50, UnspecialNum = 10):
		filter_IDs = []
    
		for i in blastScore:
			Alls = blastScore[i]
			count_S = Counter(Alls)
			x = sorted(count_S.items(), key= lambda x: x[0], reverse= True)       # shall not use count_S.most_common(), because we will sort it by score, but not by Number 
        
			c = x[1:6]        # sum all un special blasted sequences.
        
			bestHitN = x[0][1]       # number of hsps with best hit e-value
        
			bestHitScore = x[0][0]
        
			if len(c) > 1:
				testScoreDrop =  x[0][0]  - x[1][0]   # how much does second best score less than best score 
			else:
				testScoreDrop = ScoreDrop * 2
        
			SubNum = sum([i[1] for i in c if i[0] < bestHitScore - ScoreDrop])       # subordinate hsps numbers with rank 2  -  6
        
			if bestHitN <= hitNum and (testScoreDrop > ScoreDrop and SubNum < UnspecialNum):
				filter_IDs.append(i)
            
		outN = self.tmp + "_filter_ProbeHits.xls"
		inH = open(self.tmp, "r")
    
		outH = open(outN, "w")
		for line in inH:
			HFields = line.strip().split("\t")
			if HFields[0] in filter_IDs:
				outH.write(line)
            
		outH.close()
		inH.close()
    
		return filter_IDs
    
	def filter_selfDimer(self, DeltaG_dict, dG = -5.6):
		pass
		# return list
	def combine_Blastn_dG(self, blastList, dimerList, outFile):
		pass
		# return None

if __name__ == '__main__':

	import argparse
	
	HELP = """USAGE: python {0} -i sample_name """.format(__file__)
             
	parser = argparse.ArgumentParser( description = HELP ) 
 
	parser.add_argument("-d", "--dimerL", action="store", type=int, default=14, help="hetero dimer maximal length")
	parser.add_argument("-i", "--fileName", action="store", default = "test_seqfile.txt", help="analysis seq file")
	args = parser.parse_args()
	
	#maxLen, minLen = 40, 25
	
	dL = args.dimerL
	probefile = args.fileName
	
	if not os.path.exists(probefile):
		sys.exit("    Probe information file not found. \n    Please check your filename\n")
	
	SeqsDict = get_probeSeq(probefile)
	
	for sname in SeqsDict:

		PLabel = sname
		seq    = SeqsDict[sname]
		if not os.path.exists(PLabel):
			os.mkdir(PLabel)
		os.chdir(PLabel)
		SeqManipu = ProbeSolution(seq)    ## this object will be used at do_Blastn

		Probe_attFile = "Step1_" + PLabel +".xls"
		out_probeInfo = open(Probe_attFile, "w")
		InTmRange, Tm005M, Tm1M = ProbeSolution.TemMelt(seq, 53, 60)
		#if InTmRange:
		seq_len = str(len(seq))
			#lastS = seq[-16:]
		out_probeInfo.write("{0}\t{1}\t{2}\t{3:.2f}\t{4:.2f}\t{5:.2f}\t".format(PLabel, seq_len , seq, Tm005M, Tm1M, ProbeSolution.GCCalc(seq)))
			#out_test.write("{0}_{1}_{2}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\n".format(PLabel, istart, iend, seq_len, s, Tm005M, Tm1M, ProbeSolution.GCCalc(s), lastS, ProbeSolution.GCCalc(lastS)))

		
		## create input file for heterodimer
		HetDimer_srcFile = open(sname + "_dimerHet.txt", "w")
		for sn_each in SeqsDict:
			if sn_each != sname:
				HetDimer_srcFile.write(PLabel+":" + sn_each+ " " +seq + " " + SeqsDict[sn_each] + "\n")
		
		HetDimer_srcFile.close()
		subprocess.call("/home/qiaozy2/Software/C_test/NewHeteroDimer_cal {0} {1} > {2}".format(sname + "_dimerHet.txt", dL, sname+ "_heterodimer.out"), shell=True)
		
		#generate_fa_for_blast
		with open(sname + "_blast_input.fa","w") as blastqueryf:
			blastqueryf.write(">{0}\n{1}\n".format(PLabel, seq))
		
		# generate input file for self dimer prediction
		with open(sname + "_selfDimer_input.txt","w") as dimerseq:
			dimerseq.write("{0}\t{1}\n".format(PLabel, seq))
            
        # calculate selfdimer DeltaG
        
		subprocess.call("/home/qiaozy2/Software/C_test/Test_selfDimer2 {0} > {1}".format(sname + "_selfDimer_input.txt", "Step3_"+sname+"_selfDimer.out"), shell=True)
	
		# do blastn, because the parameters -K -v -b could not change the maximal HSP number of blastall results, so we use blast+ .
		SeqManipu.do_Blastn("/home/qiaozy2/Software/ncbi-blast-2.6.0+/bin/blastn", sname + "_blast_input.fa")
	
		   

		#parse blast results, use blast+ , so we discard blastall suffix name
		#ParserA = FileParser(Probe_attFile.strip(".xls") +".fulllen_blastall.xml")
		ParserA = FileParser(sname + "_blast_input.xml")
		ParserA.writeHeader("xml",  "Step4_blastScore_" + PLabel + ".xls" )
    
		ParserA.parse_xml("Step4_blastScore_" + PLabel + ".xls")

		# get blast score     
		outN, D1 = FileParser.GetBlastScore("Step4_blastScore_" + PLabel + ".xls")
		out_probeInfo.write(",".join([str("{0:.2f}".format(i)) for i in D1[sname]]))
		out_probeInfo.write("\t")
		# filter with blast
		filter_Init = Table_Filter("Step4_blastScore_" + PLabel + ".xls")
		outlist = filter_Init.filter_blast( D1 )
	
			# get hetero dimer deltaG
		resultsDeltaG = FileParser.GetDeltaG(sname+ "_heterodimer.out")
		# print(resultsDeltaG)
		sortDG = sorted(resultsDeltaG.values(), reverse=True)
		out_probeInfo.write(",".join([str(i) for i in sortDG]))
		out_probeInfo.write("\n")
		out_probeInfo.close()
		# filter with dG

		#do combine_Blastn_dG evaluation
		os.chdir("../")
