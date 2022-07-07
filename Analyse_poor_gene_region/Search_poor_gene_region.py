#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 17 09:57:51 2022

@author: christophe

Script that relies on gene localization (ensembl.org) to find regions without 
a gene.
The analysis of regions is done on areas that are at least 3Mb away from 
telomeres and centromeres (variable size), and eliminates all regions that are 
between 2 exons of the same gene.
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
        elif previous_pos != None:
            next_pos = line.split('\t')
            if previous_pos[0] == next_pos[0]:
                diff = (int(next_pos[1]) - int(previous_pos[2]))
                if diff >= 50000 and previous_pos[3] != next_pos[3]:
                    temp = {'Region':{}}
                    temp['Region']['Chr'] = previous_pos[0]
                    temp['Region']['length'] = diff
                    temp['Region']['Start'] = int(previous_pos[2]) 
                    temp['Region']['End'] = int(next_pos[1])
                    temp['Region']['previousGene'] = previous_pos[3] 
                    temp['Region']['NextGene'] = next_pos[3]
                    poor_gene_coord.append(temp)
                previous_pos = next_pos.copy()
            elif previous_pos[0] != next_pos[0] and int(line.split('\t')[1]) >= size:
                previous_pos = next_pos.copy()


physical_domain_file = "Physical_domain_droso_DM6.csv"
folder = "/home/christophe/Documents/Projets/Droso_genetic/Data_from_Sexton"
physical_domain_path = folder + os.sep + physical_domain_file 
physical_domain_list = [] 
with open(physical_domain_path, mode='r') as file:
    for line in file:
        if line.split(",")[1] != "#VALEUR !":
            temp = list(line.replace("\n","").split(","))
            temp[1] = int(temp[1])
            temp[2] = int(temp[2])
            physical_domain_list.append(temp)
domainType = [x[3] for x in physical_domain_list]
domainType = list(set(domainType))

domain = {}
for item in domainType:
    domain[item]=0    
for item in poor_gene_coord:
    item["Domain"]=domain.copy()
    
for item in poor_gene_coord:
    for domain in physical_domain_list:
        if item['Region']['Chr'] == domain[0] :
            item_start = item['Region']['Start']
            item_end = item['Region']['End']
            if domain[1] <= item_start and item_start <= domain[2]:
                if domain[2] >= item_end:
                    item['Domain'][domain[3]] = item['Region']['length']
                else:
                    item['Domain'][domain[3]] = item['Domain'][domain[3]] + (domain[2] - item_start)
            elif item_start < domain[1] and domain[1] < item_end:
                if domain[2] <= item_end:
                    item['Domain'][domain[3]] = item['Domain'][domain[3]] + (domain[2] - domain[1])
                else:
                    item['Domain'][domain[3]] = item['Domain'][domain[3]] + (item_end - domain[1])

            
            
#%% Mise en forme de fichier .csv            

resultFolder = "/home/christophe/Python_Results"
fileFolder = "Poor_gene_region"
fileName = "Liste_poorGeneRegion.csv"
finalFolder = resultFolder + os.sep + fileFolder
pathFile = finalFolder + os.sep + fileName
if not os.path.isdir(finalFolder):
    os.mkdir(finalFolder)

with open(pathFile, mode='w') as file:
    entete = "Chr.,Region length,Start,End,Gene before region, Gene after region, "
    file.write(entete + '\n')
    for region in poor_gene_coord:
        file.write(','.join(str(n) for n in region) + '\n')
        
        