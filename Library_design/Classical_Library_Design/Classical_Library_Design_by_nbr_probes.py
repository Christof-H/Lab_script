#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 16 16:07:15 2021

@author: christophe

Ce Script permet de choisir une/des régions et de faire les sondes primaires
avec 1 bcd/RT différent par région différente (avec le choix possible de bcd/RT par 
sonde primaire de 2 à 5 pour le même bcd/RT)

"""
import os, json, sys

###-------------------IMPORTATION DES PARAMETRES DE LA LIBRAIRIE------------------------

rootFolder = os.getcwd()


jsonFile = 'input_parameters.json'
jsonPath = rootFolder + os.sep + jsonFile
with open(jsonPath, mode='r') as file:
    input_parameters = json.load(file)


chromosomeFile = input_parameters["ChromosomeFile"]
chromosomeFolder = input_parameters["ChromosomeFolder"]
resolution = input_parameters["resolution"]
startLib = input_parameters["startLib"]
nbrLociTotal = input_parameters["nbrLociTotal"]
nbrProbeByLocus = input_parameters["nbrProbeByLocus"]
nbrBcd_RT_ByProbe = input_parameters["nbrBcd_RT_ByProbe"]
PrimerU = input_parameters["PrimerU"]
bcd_RT_File = input_parameters["bcd_RT_File"]





primerUnivFile = 'Primer_univ.csv'

bcd_RT_Path = rootFolder + os.sep + bcd_RT_File
primaryPath = chromosomeFolder + os.sep + chromosomeFile
primerUnivPath = rootFolder + os.sep + primerUnivFile



###-------------------CREATION DU DOSSIER RESULTAT----------------------
resultFolder = rootFolder + os.sep + 'Library_Design_Results'
if not os.path.exists(resultFolder):    
    os.mkdir(resultFolder)

###-----------VERIFICATION DE LA PRESENCE DE TOUS LES FICHIERS----------

# A faire ultérieurement............


#%%-----------------FORMATAGE DES FICHIERS EN VARIABLES-------------------
from functions import FormatFile


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
# Sort list of primary probes by their coordinates (element[0]):
listSeqGenomic.sort(key=lambda x: int(x[0]))

# Ouverture et formatage des primers universels dans la variable primerUniv :
primerUniv = dict()
replace_primer = ['\n']
split_primer = [',']
FormatFile (primerUnivPath, primerUniv, 'primer',split_list=split_primer,replace_list=replace_primer)


print('-'*70)
print('listSeqGenomic =', listSeqGenomic[0])
print('-'*70)
print('bcd_RT =', bcd_RT[:2])
print('-'*70)
print('primerUniv = ', 'primer1 =', primerUniv['primer1'])
print('-'*70)

#%%----REMPLISSAGE DES LOCUS (Primers Univ, start, end, Seq DNA genomic)-------
from functions import LocusDataClass
import random

# recherche des primers universels souhaités
primer = [primerUniv[x] for x in primerUniv.keys() if x == PrimerU]
primer = primer[0]


# instanciation de la classe LocusDataClass avec tous les loci et regroupement 
# des objets locus dans la liste total_locus
total_locus = list()
for i in range(nbrLociTotal):
  total_locus.append(LocusDataClass(locusN=i+1, chrName=chromosomeFile.split('.')[0], primers_Univ = primer, startSeq=0, endSeq=0))


# Réduction de la listSeqGenomic afin d'éviter d'irérer sur la liste entière pour chaque locus:
listSeqGenomicSorted=[]
count=1
for seq in listSeqGenomic :
    if startLib <= int(seq[0]) and count <= (nbrLociTotal * nbrProbeByLocus):
        listSeqGenomicSorted.append(seq)
        count += 1
          

# Attribution des sequences d'ADN complémentaires sondes primaires par locus en
# fonction de leur nombre par locus (nbrProbeByLocus)
mini = 0
maxi = nbrProbeByLocus
for locus in total_locus :
    temp=[]
    locus.startSeq = listSeqGenomicSorted[mini][0]
    locus.endSeq = listSeqGenomicSorted[maxi-1][1]
    for seq in listSeqGenomicSorted[mini:maxi] :
        temp.append(seq[2])
    locus.seqProbe = temp
    mini += nbrProbeByLocus
    maxi += nbrProbeByLocus


# Affichage pour exemple d'un locus :
locus = total_locus[0].__dict__
[print(x,':',locus[x]) for x in locus.keys()]

#%%COMPLETION DES SEQUENCES PRIMAIRES AVEC BARCODES ou RT, ET PRIMERS UNIVERSELS

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
print('-'*70)


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
#Cr&ation du dossier daté pour différentier les différents librairies dessinées
import datetime as dt
date_now = dt.datetime.now().strftime("%Y%m%d_%H%M")


pathResultFolder = resultFolder + os.sep + date_now
os.mkdir(pathResultFolder)


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
    file.write('Chromosome,Locus_N,Start,End,Region size,Barcode,PU.Fw,PU.Rev,Nbr_Probes\n')
    for locus in total_locus :
        size = str(int(locus.endSeq) - int(locus.startSeq))
        file.write(str(locus.chrName)+','+str(locus.locusN)+','+locus.startSeq\
+','+locus.endSeq+','+size+','+str(locus.bcdLocus)+','+locus.primers_Univ[0]+','\
+locus.primers_Univ[2]+','+str(len(locus.seqProbe))+'\n')  

# Sauvegarde des parametres ayant servis pour générer la bibliothèque sous 
# forme d'un fichier.json

parametersFilePath = pathResultFolder + os.sep + '4-OutputParameters.json'

#Récupérartion des paramètres
parameters = {}
parameters['Script_Name']='Classical_Library_Design.py'
parameters['ChromosomeFile']=chromosomeFile
parameters['ChromosomeFolder']=chromosomeFolder
parameters['resolution']=resolution
parameters['startLib']=startLib
parameters['endLib']=startLib+(resolution*nbrLociTotal)
parameters['nbrLociTotal']=nbrLociTotal
parameters['nbrProbeByLocus']=nbrProbeByLocus
parameters['nbrBcd_RT_ByProbe']=nbrBcd_RT_ByProbe
parameters['PrimerU']=PrimerU
parameters['bcd_TR_File']=bcd_RT_File
parameters['primerUnivFile']=primerUnivFile




# écriture du fichier json.
with open(parametersFilePath, mode="w") as file :
    json.dump(parameters, file, indent=4)

print('-'*40, '\n')
print(f'All files concerning your library design are saved in {pathResultFolder}/')
