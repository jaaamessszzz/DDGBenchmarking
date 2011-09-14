#!/usr/bin/python

import sys
import os

mutantsfile = os.path.join("rawdata", "mutants.txt")

def parseFloat(str):
	if str.strip() == "NaN":
		return None
	else:
		return float(str)
	
def main():
	F = open(mutantsfile)
	for line in F.read().split("\n"):
		if line.strip() and line[0] != '#':
			print(line)
			index = int(line[0:4])
			pdbID = line[5:9]
			chainID = line[9]
			wildtypeAA = line[14]
			mutantAA = line[19]
			resID = int(line[20:25])
			conformation = line[25]
			experimentalScore = parseFloat(line[26:35])
			numMeasurements = int(line[36:41])
			
			CCPBSA = (parseFloat(line[41:51]), int(line[51:53]) == 1)
			EGAD = (parseFloat(line[53:63]), int(line[63:65]) == 1)
			FoldX = (parseFloat(line[65:75]), int(line[75:77]) == 1)
			Hunter = (parseFloat(line[77:87]), int(line[87:89]) == 1)
			IMutant2 = (parseFloat(line[89:99]), int(line[99:101]) == 1)
			Rosetta = (parseFloat(line[101:111]), int(line[111:113]) == 1)
			
			print(index,pdbID,chainID,wildtypeAA,mutantAA,resID,experimentalScore,numMeasurements, CCPBSA, EGAD, FoldX, Hunter, IMutant2, Rosetta)
	F.close()
	
main() 
