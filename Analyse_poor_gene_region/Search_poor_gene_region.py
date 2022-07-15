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


# Import Physical Domain form Sexton paper

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
    
for region in poor_gene_coord:
    for domain in physical_domain_list:
        if region['Region']['Chr'] == domain[0] :
            region_start = region['Region']['Start']
            region_end = region['Region']['Start']
            if domain[1] <= region_start and region_start <= domain[2]:
                if domain[2] >= region_end:
                    region['Domain'][domain[3]] = region['Region']['length']
                else:
                    region['Domain'][domain[3]] = region['Domain'][domain[3]] + (domain[2] - region_start)
            elif region_start < domain[1] and domain[1] < region_end:
                if domain[2] <= region_end:
                    region['Domain'][domain[3]] = region['Domain'][domain[3]] + (domain[2] - domain[1])
                else:
                    region['Domain'][domain[3]] = region['Domain'][domain[3]] + (region_end - domain[1])


# Import enhancer from RedFly site

# Search for Cis-regulatory elements in identified gene-poor regions
# CMRs: Cis regulatory modules
# CMRs segment (CMRS_seg): large CMR DNA region not really defined
# pCMRs: predicted CMRs
# RCs: Reporter construct
# TFBSs: Transcription factor binding sites

for item in poor_gene_coord:
    item['CMR']={'CMRs':[], 'CMR_seg':[], 'pCMRs':[], 'RCs':[], 'TFBSs':[]}

cmr_folder = "/home/christophe/Documents/Projets/Droso_genetic/RedFly_data"
cmrs_file = "all_drosophila_melanogaster_crms.bed"
cmrsSeg_file = "all_drosophila_melanogaster_crm_segments.bed"
pCMRs_file = "all_drosophila_melanogaster_predicted_crms.bed"
rcs_file = "all_drosophila_melanogaster_rcs.bed"
tfbss_file = "all_drosophila_melanogaster_tfbss.bed"

cmrs_path = cmr_folder + os.sep + cmrs_file
cmrsSeg_path = cmr_folder + os.sep + cmrsSeg_file
pCMRs_path = cmr_folder + os.sep + pCMRs_file
rcs_path = cmr_folder + os.sep + rcs_file
tfbss_path =cmr_folder + os.sep + tfbss_file

cmrs_list = []
cmrsSeg_list = []
pCMRs_list = []
rcs_list = []
tfbss_list = [] 

def create_CMR_list(doc_path, cmr_list):
    with open(doc_path, mode='r') as file:
        for line in file:
            temp = list(line.replace("\n","").split('\t'))
            temp[1] = int(temp[1])
            temp[2] = int(temp[2])
            cmr_list.append(temp)

def search_CMR(cmr_list, cmr_name, dic):
    for region in dic:
        temp=[]
        for item in cmr_list:
            if item[0] == region['Region']['Chr'] \
                and (region['Region']['Start']<item[1]<region['Region']['End'] \
                     or region['Region']['Start']<item[2]<region['Region']['End']):
                temp.append(item[3])
        region['CMR'][cmr_name]= temp


create_CMR_list(cmrs_path, cmrs_list)
create_CMR_list(cmrsSeg_path, cmrsSeg_list)
create_CMR_list(pCMRs_path, pCMRs_list)
create_CMR_list(rcs_path, rcs_list)
create_CMR_list(tfbss_path, tfbss_list)


search_CMR(cmrs_list, "CMRs", poor_gene_coord)
search_CMR(cmrsSeg_list, "CMR_seg", poor_gene_coord)
search_CMR(pCMRs_list, "pCMRs", poor_gene_coord)
search_CMR(rcs_list, "RCs", poor_gene_coord)
search_CMR(tfbss_list, "TFBSs", poor_gene_coord)


# Import MiMIC strain coordinate

mimicFolder = "/home/christophe/Documents/Projets/Droso_genetic/MiMIC_line_from_Flybase"
mimicfile = "FlyBase_Fields_download.txt"
mimicpath = mimicFolder + os.sep + mimicfile
mimic_list = []

with open(mimicpath, mode='r') as file:
    lines = file.readlines()[1:] #remove header of text file
    for line in lines:
        temp=[]
        temp = line.replace('\n', '').split('\t')
        coord = temp[6].replace(':', ' ').replace('..', ' ').split(' ')
        if temp[6] != '-':      #remove MiMIC strain without coordinate
            mimicStrain=[]
            coord[1] = int(coord[1])
            coord[2] = int(coord[2])
            mimicStrain.extend(coord)
            mimicStrain.extend([temp[0], temp[4], temp[9], temp[10]])
            mimic_list.append(mimicStrain)

# Check for MiMIC lines in the different gene-poor regions selected
for locus in poor_gene_coord:
    mimic_line=[]
    for mimic in mimic_list:
        if locus['Region']['Chr'] == 'chr'+mimic[0]:
            if locus['Region']['Start'] <= mimic[1] and mimic[2] <= locus['Region']['End']:
                mimic_line.append(mimic)
    locus['MiMIC'] = mimic_line

                       
            
#%% Mise en forme de fichier .csv            

resultFolder = "/home/christophe/Python_Results"
fileFolder = "Poor_gene_region"
fileName = "Liste_poorGeneRegion.csv"
finalFolder = resultFolder + os.sep + fileFolder
pathFile = finalFolder + os.sep + fileName
if not os.path.isdir(finalFolder):
    os.mkdir(finalFolder)

with open(pathFile, mode='w') as file:
    entete = "Chr.,Region length,Start,End,Gene before region, \
        Gene after region, Epigenetic class, CRMs, CRM segment, \
            RCs, Predited CRM, TFBSs, Nber MiMIC, MiMIC Line ID"
    file.write(entete + '\n')
    for region in poor_gene_coord:
        chrom = region['Region']['Chr']
        length = str(region['Region']['length'])
        start = str(region['Region']['Start'])
        end = str(region['Region']['End'])
        previous = region['Region']['previousGene']
        following = region['Region']['NextGene']
        nbCMRs = str(len(region['CMR']['CMRs']))
        nbCMR_seg = str(len(region['CMR']['CMR_seg']))
        nbpCMRs = str(len(region['CMR']['pCMRs']))
        nbRCs = str(len(region['CMR']['RCs']))
        nbTFBSs = str(len(region['CMR']['TFBSs']))
        nbMimic =str(len(region['MiMIC']))
        lineID =[]
        for mimicline in region['MiMIC']:
            lineID.append(mimicline[-1])
        
        
        list_marks =[]
        for key, value in region['Domain'].items():
            if value != 0:
                list_marks.append(key)
        domain = '/'.join(list_marks)
        file.write(chrom +','+ length +','+ start +','+ end +','+ previous \
                   +','+ following +','+ domain +','+ nbCMRs +','+ nbCMR_seg \
                   +','+ nbRCs +','+ nbpCMRs +','+ nbTFBSs +','+ nbMimic \
                    +','+' - '.join(lineID) + '\n')
        
        
        
        
        
        
        
        
        
        