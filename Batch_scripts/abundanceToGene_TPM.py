#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 12 14:07:57 2019

Converting Kallisto abundances into a gene matrix
@author: alim
"""

import sys, os
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
#import seaborn as sns
#import numpy as np


if len(sys.argv) < 4:
    print('usage: %s expID geneTranscriptFile abundancesPath' % sys.argv[0])
    sys.exit(0)

expID = sys.argv[1]
geneFilename = sys.argv[2]
abundancesPath = sys.argv[3]

geneFile = open(geneFilename,'r')

geneIDList = []
cellList = []
cellIDList = []
cellSortList = []
transcriptIDTOgeneIDDict = {}

print("processing genes")
lines = geneFile.readlines()
for line in lines:
    fields = line.strip().split()
    geneENSG = fields[0]
    geneName = fields[1]
    geneID = '%s-%s' % (geneENSG,geneName)
    transcriptENST = fields[2]
    transcriptName = fields[3]
    if geneID not in geneIDList:
        geneIDList.append(geneID)
    transcriptIDTOgeneIDDict[transcriptENST] = geneID
geneFile.close()

print("processing abundances")  
for filename in os.listdir(abundancesPath):
    if filename.endswith('.tsv'):
        (cell,suffix) = filename.split('.')
        cellID = cell
        cellName = expID + '_' + str(cellID)
        cellSortList.append((cellID,cellName, abundancesPath + '/' + filename))
cellSortList.sort()

for (cellID, cellName, cellFilename) in cellSortList:
        #print(cell)
        cellList.append(cellName)
        cellIDList.append(cellID)

countDF = pd.DataFrame(index=geneIDList,columns=cellList)
countDF = countDF.fillna(0.)

cellCounter = 0
for (cellID, cellName, cellFilename) in cellSortList:
        print(cellFilename)
        abundanceFile = open(cellFilename)
        lines = abundanceFile.readlines()
        for line in lines[1:]:
            fields = line.strip().split()
            tpm = float(fields[4])
            if tpm > 0:
                transcriptID = fields[0].split('|')[0]
                geneID = transcriptIDTOgeneIDDict[transcriptID]
                countDF[cellName][geneID] += tpm
        abundanceFile.close()
        cellCounter += 1
        if cellCounter % 100 == 0:
            print(cellCounter)
print('processed %d cells' % cellCounter)

countDF = countDF.round()
countDF = countDF.astype(int)
    
print('writing csv file %s.matrix.csv.gz' % expID)
countDF.to_csv(expID + '.matrix.csv.gz',compression='gzip')
