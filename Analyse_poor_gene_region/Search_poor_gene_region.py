#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 17 09:57:51 2022

@author: christophe

Script that relies on gene localization (ensembl.org) to find regions without a gene.

"""

import os

rootFolder = "/mnt/PALM_dataserv/DATA/Commun/genomes/RefSeqs"
exonFile = "dm6.ensGene_exon.bed"

list_exon_file = rootFolder + os.sep + exonFile
#Size of the unanalyzed region extending from the beginning of the telomere or centromere
size = 3000000
poor_gene_coord = []

with open(list_exon_file, mode='r') as file:
    previous_pos = None
    for line in file:
        if previous_pos == None and int(line.split('\t')[1]) >= size:
            previous_pos = line.split('\t')
            continue
        else:
            next_pos = line.split('\t')
        if previous_pos != None:
            if previous_pos[0] == next_pos[0]:
                diff = (int(next_pos[1]) - int(previous_pos[2]))
                previous_pos = next_pos.copy()
                temp =[previous_pos[0], diff, previous_pos[2], next_pos[1], previous_pos[3]]
                if diff >= 50000:
                    poor_gene_coord.append(temp)
            elif previous_pos[0] != next_pos[0] and int(line.split('\t')[1]) >= size:
                previous_pos = next_pos.copy()



resultFolder = "/home/christophe/Python_Results"
fileFolder = "Poor_gene_region"
fileName = "Liste_poorGeneRegion.csv"
finalFolder = resultFolder + os.sep + fileFolder
pathFile = finalFolder + os.sep + fileName
if not os.path.isdir(finalFolder):
    os.mkdir(finalFolder)

with open(pathFile, mode='w') as file:
    entete = "Chr.,Region length,Start,End,Gene"
    file.write(entete + '\n')
    for region in poor_gene_coord:
        file.write(','.join(str(n) for n in region) + '\n')
        
        