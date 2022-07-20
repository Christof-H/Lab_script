#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 16 16:07:15 2021

@author: christophe

Ce Script permet de choisir une/des régions et de faire les sondes primaires
avec 1 bcd/RT différent par région différente (avec le choix possible de bcd/RT par 
sonde primaire de 2 à 5 pour le même bcd/RT)

"""
import os, glob

###-------------------PARAMETRES DE LA LIBRAIRIE------------------------

chromosome = 'chr2L_original' #'chr2L' or 'chr2R' or 'chr3L' or 'chrX' or 'ChrY'.....
resolution = 3000 # Taille des loci en nucléotides
startLib = 16411158 # Coordonnée génomique du début du 1er locus
nbrProbeByLocus = 45 # nombre de sondes primaires par locus
nbrLociTotal = 22 # nombre de loci au total
PrimerU = 'primer1' # choix du couple de primers universels 'primer1', 'primer2' ou 'primer3'....
nbrBcd_RT_ByProbe = 3 # Nbre de même barcode par sonde primaire max=5

###--------CREATION DES CHEMINS D'ACCES POUR LES FICHIERS-----------
#chemin d'accès pour le dossier Combinatorial_Library_Design
folderChromosome = '/mnt/PALM_dataserv/DATA/Commun/genomes/dm3/OligoMiner/dm3_32nt'
rootFolder = '/home/christophe/Documents/Repositories/Lab_script/Library_design'
chromosomeFile = chromosome + '.bed'
bcd_RT_File = 'List_RT.csv' # Choix entre 'Barcodes.csv' ou 'List_RT.csv'
primerUnivFile = 'Primer_univ.csv'

bcd_RT_Path = rootFolder + os.sep + bcd_RT_File
primaryPath = folderChromosome + os.sep + chromosomeFile
primerUnivPath = rootFolder + os.sep + primerUnivFile

os.chdir(rootFolder)

###-------------------CREATION DU DOSSIER RESULTAT----------------------
resultFolder = os.path.expanduser('~/Python_Results')
Lib_Design_Folder = 'Library_Design_Diff'
pathResultFolder= resultFolder+ os.sep + Lib_Design_Folder
if not os.path.exists(resultFolder):
    os.mkdir(resultFolder)
if not os.path.exists(pathResultFolder):    
    os.mkdir(pathResultFolder)

###-----------VERIFICATION DE LA PRESENCE DE TOUS LES FICHIERS----------

# A faire ultérieurement............


#%%-----------------FORMATAGE DES FICHIERS EN VARIABLES-------------------
from functions import FormatFile
import sys

# Ouverture et formatage des barcodes dans la variable barcodes :
replace_bcd_RT =['\n']
split_bcd_RT = [',']
bcd_RT= list()
FormatFile (bcd_RT_Path, bcd_RT, 'bcd_RT',split_list=split_bcd_RT, replace_list=replace_bcd_RT)

# Ouverture et formatage des coordonnées et des séquences des sondes primaires 
# dans la variable listSeqGenomic :
listSeqGenomic = list()
split_primary = ['\t']
FormatFile (primaryPath, listSeqGenomic, 'primaryProbe',split_list=split_primary)

# Ouverture et formatage des primers universels dans la variable primerUniv :
primerUniv = dict()
replace_primer = ['\n']
split_primer = [',']
FormatFile (primerUnivPath, primerUniv, 'primer',split_list=split_primer,replace_list=replace_primer)


print('-'*70)
print('listSeqGenomic =', listSeqGenomic[0])
print('-'*70)
print('barcodes =', bcd_RT[:2])
print('-'*70)
print('primerUniv = ', 'primer1 =', primerUniv['primer1'])

#%%----REMPLISSAGE DES LOCUS (Primers Univ, start, end, Seq DNA genomic)-------
from functions import LocusDataClass
import random


# recherche des primers universels souhaités
primer = [primerUniv[x] for x in primerUniv.keys() if x == PrimerU]
primer = primer[0]

# Calcul des start et end position de chaque locus :
startPositions = [startLib + x*resolution for x in range(nbrLociTotal)]
endPositions = [startLib + (x+1)*resolution for x in range(nbrLociTotal)]

# Remplissage de la classe LocusDataClass avec tous les Loci nécécesaires
total_locus = list()
for i, start, end in zip(range(nbrLociTotal),startPositions,endPositions):
  total_locus.append(LocusDataClass(locusN=i+1, chrName=chromosome, startSeq=start, endSeq=end, primers_Univ = primer))


# Attribution des sequences d'ADN complémentaires sondes primaires par locus en fonction des coordonnées génomiques
for locus in total_locus :
    temp = []
    for seq in listSeqGenomic :
            if locus.startSeq <= int(seq[0]) and int(seq[1])< locus.endSeq:
                temp.append([seq[0],seq[2]])                    
            else :
                pass
    
    if len(temp) > nbrProbeByLocus:
        random.shuffle(temp)
        temp = temp[:nbrProbeByLocus]
        temp.sort()
    locus.seqProbe = [x[1] for x in temp]

# Affichage pour exemple d'un locus :
locus = total_locus[0].__dict__
[print(x,':',locus[x]) for x in locus.keys()]

#%%COMPLETION DES SEQUENCES PRIMAIRES AVEC BARCODES ET PRIMERS UNIVERSELS

import copy
# Insertion des binding sites pour les barcodes des loci selon le schéma suivant :
# Bcdx_SeqADNgenomic_Bcdx_ ou Bcdx_SeqADNgenomic_Bcdx_Bcdx
count = 0
for locus in total_locus :
    locus.bcdLocus = bcd_RT[count][0]
    seqWithBcd = []
    if nbrBcd_RT_ByProbe == 2 :
        for item in locus.seqProbe :
            seqWithBcd.append(bcd_RT[count][1]+' '+ item +' '+bcd_RT[count][1])
        count +=1
        locus.seqProbe = seqWithBcd
    elif nbrBcd_RT_ByProbe == 3 :
        for item in locus.seqProbe :
            seqWithBcd.append(bcd_RT[count][1]+' '+ item +' '+bcd_RT[count][1]*2)
        count +=1
        locus.seqProbe = seqWithBcd
    elif nbrBcd_RT_ByProbe == 4 :
        for item in locus.seqProbe :
            seqWithBcd.append(bcd_RT[count][1]*2+' '+ item +' '+bcd_RT[count][1]*2)
        count +=1
        locus.seqProbe = seqWithBcd
    elif nbrBcd_RT_ByProbe == 5 :
        for item in locus.seqProbe :
            seqWithBcd.append(bcd_RT[count][1]*3+' '+ item +' '+bcd_RT[count][1]*2)
        count +=1
        locus.seqProbe = seqWithBcd
    

# Insertion des primers universels selon le schéma suivant:
# pu.fw_(Bcd-region_Bcdx_Bcdy_SeqADNgenomic_Bcdx_Bcdy)_pu.rev

for locus in total_locus :
    pFw=copy.deepcopy(locus.primers_Univ[1])
    pRev=copy.deepcopy(locus.primers_Univ[3])
    temp = list()
    temp=[pFw+' '+ x +' '+pRev for x in locus.seqProbe]
    locus.seqProbe = temp

# Affichage pour exemple d'une séquence de sonde primaire :    
print('-'*70)
print("exemple d'une séquence primaire :")
print('-'*70)
print(total_locus[0].seqProbe[0])


#%%------------VERIFICATION TAILLE DES SONDES PRIMAIRES-------------------
#--------------ET COMPLETION SI TAILLE TROP DIFFERENTE--------------------
from functions import Check_Length_Seq_Diff
from functions import Completion

# Cacul de la différence de taille entre les séquences des sondes primaires
diff_pourcent, max_seq_length = Check_Length_Seq_Diff(total_locus)

# Completion des séquences avec nucléotides aléatoires
# ATTENTION: complétion en 3' de la séquence !!!!!
Completion(diff_pourcent,max_seq_length,total_locus)
print('-'*70)
print("exemple de séquences primaires :")
print(total_locus[0].seqProbe[:3])

#%%----------ECRITURE DES DIFFERENTS FICHIERS RESULTATS---------------------            

#fichier détaillé avec information et séquences : Library_details
resultDetails = pathResultFolder+os.sep+'1_Library_details.txt'
with open (resultDetails, 'w') as file :
    for locus in total_locus :
        file.write('Chromosome:'+str(locus.chrName)+' Locus_N°'+str(locus.locusN)\
+' Start:'+str(locus.startSeq)+' End:'+str(locus.endSeq)+' Bcd_locus:'+locus.bcdLocus+'\n')
        for seq in locus.seqProbe :
            file.write(seq+'\n')
            
#fichier avec toutes les séquences (sans espace) uniquement : Full_sequence_Only
fullSequence = pathResultFolder+os.sep+'2_Full_sequence_Only.txt'
with open (fullSequence, 'w') as file :
    for locus in total_locus :
        for seq in locus.seqProbe :
            file.write(seq.replace(' ','')+'\n')
            
#fichier avec résumé des informations (sans séquence) : Library_Summary
Summary = pathResultFolder+os.sep+'3_Library_Summary.csv'
with open (Summary, 'w') as file :
    file.write('Chromosome,Locus_N°,Start,End,Barcode,PU.Fw,PU.Rev,Nbr_Probes\n')
    for locus in total_locus :
        file.write(str(locus.chrName)+','+str(locus.locusN)+','+str(locus.startSeq)\
+','+str(locus.endSeq)+','+str(locus.bcdLocus)+','+locus.primers_Univ[0]+','\
+locus.primers_Univ[2]+','+str(len(locus.seqProbe))+'\n')  

# Sauvegarde des parametres ayant servis pour générer la bibliothèque sous 
# forme d'un fichier.json
from functions import SaveJson

parameters = {}
parameters['Script_Name']='Classical_Library_Design.py'
parameters['chromosomeFile']=chromosomeFile
parameters['resolution']=resolution
parameters['startLib']=startLib
parameters['endLib']=startLib+(resolution*nbrLociTotal)
parameters['nbrLociTotal']=nbrLociTotal
parameters['nbrProbeByLocus']=nbrProbeByLocus
parameters['nbrBcd_RT_ByProbe']=nbrBcd_RT_ByProbe
parameters['PrimerU']=PrimerU
parameters['bcd_TR_File']=bcd_RT_File
parameters['primerUnivFile']=primerUnivFile

parametersFilePath = pathResultFolder + os.sep + '4-OutputParameters.json'

SaveJson(parametersFilePath,parameters)



