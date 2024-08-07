#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 31 09:21:04 2021

@author: christophe

Script qui permet le dessin d'une librairie avec des barcodes et/ou des MER avec
ou sans fiducial (MER) avec possibilité de nombre de sondes primaires diiférent
Exemple : 
pu.fw-MER_Fiducial-MER-Bcd-SeqADNg-MER-Bcd-pu.rev
ou
pu.fw-MER_Fiducial-Bcd-SeqADNg-Bcd-pu.rev
ou
pu.fw-MER_Fiducial-MER-SeqADNg-MER-pu.rev
ou
pu.fw-Bcd-SeqADNg-Bcd-pu.rev
ou
pu.fw-MER-SeqADNg-MER-pu.rev
        
"""

import os, sys

###-------------------PARAMETRES DE LA LIBRAIRIE------------------------

chromosome = 'chrX' #'chr2L' or 'chr2R' or 'chr3L' or 'chrX'.....
resolution = 40000 # Taille des loci en nucléotides
startLib = 9000000 # Coordonnée génomique du début du 1er locus

# Nombre de sondes primaires par locus, si vous
# voulez le même nombre de sondes primaires pour tous les loci ex:300 alors rentrez
# nbrProbeByLocus = [300]
nbrProbeByLocus = [90, 150, 300, 400] 

nbrLociTotal = 20 # nombre de loci au total
PrimerU = 'primer1' # choix du couple de primers universels 'primer1', 'primer2' ou 'primer3'....
nbrBS_ByProbe = 5 # Nbre de binding site par sonde primaire pour MER ou barcodes

fiducial= True # Mettre False si vous ne voulez pas de fiducial
list_fiducial = ['mer29', 'gttgattcaaagatccggcg']

# Indiquer ici si vous voulez utiliser que les MER ou les barcodes dans les sondes primaires
MER_or_bcd = 'MER' # entrez 'bcd' pour n'utiliser que des barcodes, si vous voulez utilisez
# les MER + bcd alors peu importe (cela sera déterminé avec use_of_bcd_with_MER)


# Indiquer ici si vous voulez utiliser les barcodes en même temps que les MER car
# il faut supprimer les barcodes qui ont des séquences similiaires avec les MER:
use_of_bcd_with_MER = True  # Mettre False si vous ne voulez pas combiner barcode + MER dans le dessin  
# MER utilisé comme imaging oligo avec les barcodes et qui ne doit pas être présent
# si on utilise un marquage combiné de MER et barcode
imaging_oligo_Bcd = ['mer1', 'caccgacgtcgcatagaacg']



###--------CREATION DES CHEMINS D'ACCES POUR LES FICHIERS-----------
#chemin d'accès pour le dossier Library_Design_DiffNbrProbeByLoci_MER_Bcd
folderChromosome = '/mnt/PALM_dataserv/DATA/Commun/genomes/dm6/OligoMiner/dm6_balanced'
# Le rootFolder est le dossier contenant les fichiers functions, barcodes....
rootFolder = '/home/christophe/Documents/Informatique/Python/Scripts/Python_Script_Tof_git/Library_design'
chromosomeFile = chromosome + '.bed'
barcodeFile = 'Barcodes.csv'
primerUnivFile = 'Primer_univ.csv'
merFile = 'seq_MER_on_primary.csv'
barcodePath = rootFolder + os.sep + barcodeFile
merPath = rootFolder + os.sep + merFile
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

# Permet de répéter les nombres de sondes primaires si différents d'un locus a l'autre
# et de vérifier qu'il y a un multiple entier entre nbrProbeByLocus et nbrLociTotal
if nbrLociTotal % len(nbrProbeByLocus) != 0 :
    sys.exit('Il ne peut y avoir un multiple de nbrProbeByLocus pour un nombre total de nbrLociTotal')
nbrProbeByLocusFinal =nbrProbeByLocus * int((nbrLociTotal/len(nbrProbeByLocus)))

# Ouverture et formatage des barcodes dans la variable barcodes :
replace_barcode =['\n']
split_barcode = [',']
barcodes= list()
FormatFile (barcodePath, barcodes, 'barcode',split_list=split_barcode, replace_list=replace_barcode)
#ATTENTION si vous voulez utiliser des Barcodes et des MER ensemble, il faut supprimer
# les 20 premiers barcodes qui présentent des séquences communes avec les MER :
if use_of_bcd_with_MER :
    del barcodes[:20]


# Ouverture et formatage des RTs MER dans la variable MER :
list_MER = list()
replace_mer =['\n']
split_mer = [',']
FormatFile (merPath, list_MER, 'barcode',split_list=split_mer, replace_list=replace_mer)
#ATTENTION si vous utilisez un fiducial celui-ci doit être retirer de la liste des MER
if fiducial :
    list_MER.remove(list_fiducial)
    
    
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
print('barcodes =', barcodes[:2])
print('-'*70)
print('list_MER =', list_MER[:2])
print('-'*70)
print('primerUniv = ', 'primer1 =', primerUniv['primer1'])


#%%----REMPLISSAGE DES LOCUS (Primers Univ, start, end, Seq DNA genomic)-------
from functions import LocusDataClass
import random

# Si utilisation des MER et des barcodes dans même bibliothèque, partie qui enleve
# le MER Imaging oligo de la list_MER qui sera utiliser avec détecter les barcodes
if use_of_bcd_with_MER :
    list_MER.remove(imaging_oligo_Bcd)

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
count_locus = 0
for locus in total_locus :
    temp = []
    for seq in listSeqGenomic :
            if locus.startSeq <= int(seq[0]) and int(seq[1])< locus.endSeq:
                temp.append([seq[0],seq[2]])                    
            else :
                pass
    nbr_probe = nbrProbeByLocusFinal[count_locus]
    random.shuffle(temp)
    temp = temp[:nbr_probe]
    temp.sort()
    locus.seqProbe = [x[1] for x in temp]
    count_locus += 1

# Affichage pour exemple d'un locus :
locus = total_locus[0].__dict__
[print(x,':',locus[x]) for x in locus.keys()]

#%%COMPLETION DES SEQUENCES PRIMAIRES AVEC CODES LOCUS,CODES REGIONS ET PRIMERS UNIVERSELS

# les sondes primaires seront composées comme ci-dessous en fonction du nombre
# de binding sites (BS) pour les barcodes et/ou MER:
# PU-BS1-ADN genomic-BS2-PU
# PU-BS1-BS2-ADN genomic-BS3-PU
# PU-BS1-BS2-ADN genomic-BS3-BS4-PU
# PU-BS1-BS2-BS3-ADN genomic-BS4-BS5-PU
# si fiducial, le fiducial sera forcement en BS1



import copy
# Insertion des binding sites pour le fiducial, et/ou barcode et/ou MER des loci :
count = 0
for locus in total_locus :
    seq5prime_end = []
    seq3prime_end = []
    seqWithBcd = []
    temp_MER_Bcd = []
    if fiducial :
        seq5prime_end.append(list_fiducial[1])
        temp_MER_Bcd.append(list_fiducial[0])
        
        
    if nbrBS_ByProbe == 2 :
        if MER_or_bcd == 'MER':
            seq5prime_end.append(list_MER[count][1])
            temp_MER_Bcd.append(list_MER[count][0])
        elif MER_or_bcd == 'bcd' :
            seq5prime_end.append(barcodes[count][1])
            temp_MER_Bcd.append(barcodes[count][0])
        for item in locus.seqProbe :
            seqWithBcd.append(temp_MER_Bcd[0]+' '+ item +' '+temp_MER_Bcd[0])
        count +=1
        locus.seqProbe = seqWithBcd
        locus.bcdLocus = temp_MER_Bcd
    if nbrBS_ByProbe == 5 or nbrBS_ByProbe == 4:
        if use_of_bcd_with_MER :
            seq5prime_end.append(list_MER[count][1])
            seq5prime_end.append(barcodes[count][1])
            seq3prime_end.append(list_MER[count][1])
            seq3prime_end.append(barcodes[count][1])
            temp_MER_Bcd.append(list_MER[count][0])
            temp_MER_Bcd.append(barcodes[count][0])
        else :
            pass # definir bloc si autres besoins....

        for item in locus.seqProbe :
            seqWithBcd.append(' '.join(seq5prime_end)+' '+ item +' '+' '.join(seq3prime_end))
        count +=1
        locus.seqProbe = seqWithBcd
        locus.bcdLocus = ' '.join(temp_MER_Bcd)
    

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
resultDetails = pathResultFolder+os.sep+'1_Library_details'
with open (resultDetails, 'w') as file :
    for locus in total_locus :
        file.write('Chromosome:'+str(locus.chrName)+' Locus_N°'+str(locus.locusN)\
+' Start:'+str(locus.startSeq)+' End:'+str(locus.endSeq)+' Bcd_locus:'+locus.bcdLocus+'\n')
        for seq in locus.seqProbe :
            file.write(seq+'\n')
            
#fichier avec toutes les séquences (sans espace) uniquement : Full_sequence_Only
fullSequence = pathResultFolder+os.sep+'2_Full_sequence_Only'
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
parameters['Script_Name']='Library_Design_DiffNbrProbeByLoci_MER_Bcd.py'
parameters['chromosomeFile']=chromosomeFile
parameters['resolution']=resolution
parameters['startLib']=startLib
parameters['endLib']=startLib+(resolution*nbrLociTotal)
parameters['nbrLociTotal']=nbrLociTotal
parameters['nbrProbeByLocus']=nbrProbeByLocus
parameters['nbrBS_ByProbe']=nbrBS_ByProbe
parameters['PrimerU']=PrimerU
parameters['fiducial']=fiducial
parameters['imaging_oligo_Bcd']=imaging_oligo_Bcd 
parameters['list_fiducial']=list_fiducial
parameters['MER_or_bcd']=MER_or_bcd
parameters['use_of_bcd_with_MER']=use_of_bcd_with_MER
parameters['barcodeFile']=barcodeFile
parameters['merFile'] = merFile
parameters['primerUnivFile']=primerUnivFile

parametersFilePath = pathResultFolder + os.sep + '4-OutputParameters.json'

SaveJson(parametersFilePath,parameters)

print('-'*70)
print(f"Les fichiers résultats ont été saucegardé dans le dossier {resultFolder}")
