import sys, os, re
import pandas as pd
from Bio import SeqIO


blastx_output = "/Users/dmorales/Desktop/Beta_transrate_filtering/chimera_transrate_blastx/chimera_transrate.fa.blastx"
df = pd.read_table(blastx_output, header=None)
blastx_columns = "qseqid qlen sseqid slen frames pident nident length mismatch gapopen qstart qend sstart send evalue bitscore".strip().split(' ')
df.columns = blastx_columns

PIDENT_CUTOFF = 30 #only look at HSPs >= this percentage similarity
LENGTH_CUTOFF = 100	#only look at HSPs >= this length

def qcov(df):
	return abs(df["qend"]-df["qstart"] + 1)

#df_filtered_1 = df[(df["pident"] >= PIDENT_CUTOFF)]
#df_filtered_2 = df_filtered_1[(qcov(df_filtered_1) >= LENGTH_CUTOFF)]
hsp= df[(df["pident"] >= PIDENT_CUTOFF) & (qcov(df) >= LENGTH_CUTOFF)]
multi_hit = hsp[hsp.duplicated(["qseqid"], keep=False)]


def separated(hsp1,hsp2):
	length1 = qcov(hsp1)
	length2 = qcov(hsp2)
	start = min(hsp1["qstart"],hsp1["qend"],hsp2["qstart"],hsp2["qend"])
	end = max(hsp1["qstart"],hsp1["qend"],hsp2["qstart"],hsp2["qend"])
	overlap = length1+length2 - (end-start) + 1
	
	#value of overlap can < 0 but only the upper limit maters
	if overlap < min(60, 0.2*min(length1,length2)):
		return True
	else: return False

separated(multi_hit.iloc[0,:],multi_hit.iloc[12589,:])

def expand_range(hsp1,hsp2):
	if hsp1 == []: return hsp2
	if hsp2 == []: return hsp1
	start1,end1,start2,end2 = hsp1["qstart"],hsp1["qend"],hsp2["qstart"],hsp2["11"]
	if start1 < end1 and start2 < end2:#both forward
		start,end = min(start1,start2),max(end1,end2)
	elif start1 > end1 and start2 > end2:#both reverse
		start,end = max(start1,start2),min(end1,end2)
	#no change if new hsp is of opposite direction
	hsp1["qstart"],hsp1["qend"] = start,end
	return hsp1


