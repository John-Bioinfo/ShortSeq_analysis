#!/usr/bin/env python
# -*- coding:UTF-8 -*-

__author__ = "Zongyun Qiao"
__copyright__ = "Copyright 2017, A Biotech"
__credits__ = [
    "Zongyun Qiao"]  # remember to add yourself
__license__ = "GPL"
__version__ = "0.1-dev, 20171130"
__maintainer__ = "Zongyun Qiao"
__email__ = "gulile@yeah.net"

import math
from collections import defaultdict
from collections import Counter

def calc_Pre_score(sL, e5, mismatch, un3):
    PTs = {}
    startS = 4.5

    for i in range(1, 17):
        PTs[i] = startS 
        startS = startS / 1.09 

    for j in range(17, 50):
        PTs[j] = 1

    s = 0
    if mismatch.find("=") != -1:
        p = []
    else:        
        misPos = mismatch.split(",")
        p = [int(i) for i in misPos]

    end5 = sL - e5 +1

    for u in range(un3 +1, end5+1):
        if u not in p:
            s += PTs[u] * 2
        else:
            s -= PTs[u] * 3
            
        
    penalty_3 = 6.0
    for v in range(1, un3+1):
        s -= penalty_3
        penalty_3 = penalty_3 / 1.071
        
    return s

def calc_PscoreNP(sL, e5, mismatch, un3, startS = 3.5, S5 = 0.5, reward = 1):
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


def filter_Byscore(inputFile, blastS = 16.0):
    
    outName = inputFile + ".P1score.xls"

    outFile = open(outName , "w")
    
    scoreD = defaultdict(list)

    oF = open(inputFile, "r")

    for line in oF:
        if not line.startswith("QuerySeq_"):
            xline   = line.rstrip().split("\t")
            
            seqId = xline[0]
            
            seqLen  = int(xline[1])
            end5    = int(xline[4])
            
            un3end = int(xline[9])
            
            s = calc_PscoreNP(seqLen, end5, xline[7], un3end)
            qe = 0.6293523487909248 * seqLen * 3095693983 * math.pow(10, -0.5967479540423579 * s)
            
            if s > blastS:
            
                outFile.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\n".format(xline[0], xline[1], end5, xline[6], xline[7], xline[9],xline[15], xline[16], xline[17], xline[3], s, qe))
            
                scoreD[seqId].append(s)

    outFile.close()
    oF.close()

    return outName, scoreD

    
def filter_blast( tableFile, blastScore, hitNum = 3, ScoreDrop = 30, UnspecialNum = 2):
    filter_IDs = set()
    
    for i in blastScore:
        Alls = blastScore[i]
        count_S = Counter(Alls)
        x = sorted(count_S.items(), key= lambda x: x[0], reverse= True)       # shall not use count_S.most_common(), because we will sort it by score, but not by Number 
        
        print(x)
        c = x[1:6]        # sum all un special blasted sequences.
        
        bestHitN = x[0][1]       # number of hsps with best hit e-value
        
        bestHitScore = x[0][0]
        
        if len(c) > 1:
            testScoreDrop =  x[0][0]  - x[1][0]   # how much does second best score less than best score 
        else:
            testScoreDrop = ScoreDrop * 2
        
        SubNum = sum([i[1] for i in c if i[0] > bestHitScore - ScoreDrop])       # subordinate hsps numbers with rank 2  -  6
        
        if bestHitN <= hitNum and (testScoreDrop > ScoreDrop or SubNum < UnspecialNum):
            filter_IDs.add(i)
            
    outN = tableFile + "_filter_ProbeHits.xls"
    inH = open(tableFile, "r")
    
    outH = open(outN, "w")
    for line in inH:
        HFields = line.strip().split("\t")
        if HFields[0] in filter_IDs:
            outH.write(line)
            
    outH.close()
    inH.close()
    
    return len(filter_IDs)

filename = "Class_parser_testEGFR.xls"

outN, D1 = filter_Byscore(filename)

for i in D1:
    pass
    #print("{0} => {1}".format(i, "\t".join([str(j) for j in D1[i]])))

candi_IDNum = filter_blast(filename, D1)  

print(candi_IDNum)    
      

def calS_test(sL, e5, mismatch, un3):
    PTs = {}
    startS = 30

    for i in range(1, 21):
        PTs[i] = startS * 2.5 / 10
        startS -= 1

    startS_ext = 20

    for j in range(21, 50):
        if startS_ext > 0:
            PTs[j] = startS_ext * 1.3 / 10
            startS_ext -= 1
        else:
            PTs[j] = startS_ext * 1.3 / 10
        
    s = 0
    if mismatch.find(",") != -1:
        misPos = xline[7].split(",")
        p = [int(i) for i in misPos]
    else:
        p = []

    end5 = sL - e5

    for u in range(un3 +1, end5+1):
        if u not in p:
            s += PTs[u]
            
    return s



"""

print(calc_PscoreNP(40,1,"=",0, 3.4))

out16_All = open("O_reward1_16_All.xls","w")

out16_All.write("SID_16\tM_pattern\tRight_mis\tScore\tPenalty\n")

for a in range(30,41):

    ta = a/10.0
    outF1 = open("test_reward1_"+ str(a) + ".xls", "w")
    
    outF1.write("SID\t" + "#" * 16 + "\t####\n")   
    seqNumA = 0

    for i in range(1, 17):
        for j in range(0, 17-i):
            #print (" " * (i-1) + "|" * (17-i-j) , end = "=>")
            seqNumA += 1
            outF1.write("PA_" + str(seqNumA) + "\t" + " " * (i-1) + "|" * (17-i-j) + "\t")
            sj = calc_PscoreNP(16, i, "=", j, startS = ta)
            #print("{0}\t{1:.2f}".format(j, sj))
            outF1.write("{0}\t{1:.2f}\n".format(j, sj))
            
            out16_All.write("PA_" + str(seqNumA) + "\t" + " " * (i-1) + "|" * (17-i-j) + "\t")
            out16_All.write("{0}\t{1:.5f}\t{2}\n".format(j, sj, ta))
    outF1.close()

out16_All.close()    
    


O40_All = open("O_reward1_40_All.xls", "w")

O40_All.write("SID_40\tM_pattern\tRight_mis\tScore\tPenalty\n")

for b in range(30,41):
    tb = b/10.0
    outF2 = open("test_reward1_40L_" + str(b) + ".xls", "w")   
    outF2.write("SID\t" + "#" * 40 + "\t####\n")    
    seqNumB = 0
    
    for i in range(1, 41):
        for j in range(0, 41-i):
            #print (" " * (i-1) + "|" * (41-i-j) , end = "=>")
            
            seqNumB += 1
            outF2.write("PB_" + str(seqNumB) + "\t" + " " * (i-1) + "|" * (41-i-j)  + "\t")
            sj = calc_PscoreNP(40, i, "=", j, startS = tb)
            if i > 25:
                s_16 = calc_PscoreNP(40, i, "=", j, startS = tb)
            else:
                s_16 = calc_PscoreNP(40, 25, "=", j, startS = tb)
            #print("{0}\t{1:.2f}".format(j, sj))
            outF2.write("{0}\t{1:.2f}\t{2:.2f}\n".format(j, sj, s_16))
            
            O40_All.write("PB_" + str(seqNumB) + "\t" + " " * (i-1) + "|" * (41-i-j)  + "\t")
            O40_All.write("{0}\t{1:.2f}\t{2}\n".format(j, sj, tb))
    outF2.close()

O40_All.close()

"""

#dp1 = calc_Pscore(18, 6,"", 4)

#print ("{0}  {1}  {2}  {3}=>  {4}".format(18, 6,"", 4, dp1))

#dp2 = calc_Pscore(36, 12,"", 8)

#print ("{0}  {1}  {2}  {3}=>  {4}".format(36, 12,"", 8, dp2))
