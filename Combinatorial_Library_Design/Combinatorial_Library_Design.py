#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 16 16:07:15 2021

@author: christophe
"""
import os, glob

###-------------------PARAMETRES DE LA LIBRAIRIE------------------------

chromosome = 'chr3L' #'chr2L' or 'chr2R' or 'chr3L'.....
resolution = 40000 # Taille des loci en nucléotides
startLib = 25000 # Coordonnée génomique du début du 1er locus
nbrProbeByLocus = 90 # nombre de sondes primaires par locus
nbrRegion = 2 # nombre de régions au total
nbrSousRegion = 3
nbrLociByRegion = 30 # nombre de loci par région
nbrLoci = nbrLociByRegion*nbrRegion # nombre de loci au total (=nbre de loci par région x nbre de régions)
PrimerU = 'primer1' # choix du couple de primers universels 'primer1', 'primer2' ou 'primer3'....

###--------CREATION DES CHEMINS D'ACCES POUR LES FICHIERS-----------
#chemin d'accès pour le dossier Combinatorial_Library_Design
rootFolder = '/home/christophe/Documents/Informatique/Python/Scripts/En cours/Combinatorial_Library_Design'
codeBook_LocusFile = 'Best_Seq_Codes_30Bit_2HWeiht_51Cod_Esp10.csv'
codeBook_RegionFile = 'Best_Seq_Codes_15Bit_2HWeiht_10Cod_Esp4.csv'
chromosomeFile = chromosome + '.bed'
barcodeFile = 'Barcodes.csv'
primerUnivFile = 'Primer_univ.csv'

codeLocusPath = rootFolder + os.sep + codeBook_LocusFile
codeSousRegionPath = rootFolder + os.sep + codeBook_RegionFile
barcodePath = rootFolder + os.sep + barcodeFile
primaryPath = rootFolder + os.sep + chromosomeFile
primerUnivPath = rootFolder + os.sep + primerUnivFile

os.chdir(rootFolder)

###-----------VERIFICATION DE LA PRESENCE DE TOUS LES FICHIERS----------

# A faire ultérieurement............


#%%-----------------FORMATAGE DES FICHIERS EN VARIABLES-------------------
from functions import FormatFile
import sys

# Vérification que le nombre de loci par région est compatible avec le nombre
# de sous-region (=nombre entier) : ex : il ne peut y avoir 3 sous-région pour un 
# total de 50 loci par région (16.6 loci par sous-région)
if nbrLociByRegion % nbrSousRegion == 0 :
    pass
else :
    sys.exit('Impossible de continuer, changer le nombre de loci par région pour\
 avoir un nombre entier de loci par sous-région')    


# Ouverture et formatage des codebook contenant les codes pour les Locus 
# dans la variable codeBookLocus :
replace_list_codeBook = ["'","\n",","]    
codeBookLocus = list()
FormatFile (codeLocusPath, codeBookLocus, 'codeBook',replace_list=replace_list_codeBook)

# Ouverture et formatage des codebook contenant les codes pour les Regions 
# dans la variable codeBook_SousRegion :
replace_list_codeBook = ["'","\n",","]
codeBook_SousRegion = list()
FormatFile (codeSousRegionPath, codeBook_SousRegion, 'codeBook',replace_list=replace_list_codeBook)

# Ouverture et formatage des barcodes dans la variable barcodes :
replace_barcode =['\n']
split_barcode = [',']
barcodes= list()
FormatFile (barcodePath, barcodes, 'barcode',split_list=split_barcode, replace_list=replace_barcode)

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
print('codeBookLocus =', codeBookLocus[:2])
print('-'*70)
print('codeBook_SousRegion =', codeBook_SousRegion[:3])
print('-'*70)
print('primerUniv = ', 'primer1 =', primerUniv['primer1'])

#%%-------------------REMPLISSAGE DES LOCUS---------------------------
from functions import LocusDataClass
import random

# Calcul des start et end position de chaque locus :
startPositions = [startLib + x*resolution for x in range(nbrLoci)]
endPositions = [startLib + (x+1)*resolution for x in range(nbrLoci)]

#ATTENTION pour le remplissage des codes pour chaque locus
#la séquence  de codes pour les loci sont identique d'une région à l'autre = 
# chaque région contient une séquence de code qui est répétée dans les autres régions
#ex Région 1 : code1,code2....code30 Région 2 : code1,code2....code30 Région 3 :même chose
concatenate_CodeLocus = codeBookLocus[:nbrLociByRegion]*nbrRegion

#ATTENTION pour le remplissage des codes pour les sous-régions
#Chaque région est divisé en un certains nombre de sous-région(nbrSousRegion = 3 par défaut)
#si il y a 30 loci par région, avec 3 sous-région, il y aura donc 10 loci par
#sous-région. Ces 10 loci auront le même code région
nbrTotalSousRegion = nbrRegion * nbrSousRegion
ListFinale_codeBook_SousRegion = codeBook_SousRegion[:nbrTotalSousRegion]
concatenate_CodeSousRegion=[]
temp_CodeSousRegion =list()
for x in ListFinale_codeBook_SousRegion :
    temp_CodeSousRegion.extend([[x]*int((nbrLociByRegion/nbrSousRegion))])
for x in temp_CodeSousRegion :
    concatenate_CodeSousRegion.extend(x)


# Calcul des numéros de régions :
numeroRegion = list()
numeroRegion = [x for x in range(1,nbrRegion+1)]*nbrLociByRegion
numeroRegion.sort()

# Calcul des numéros de sous-régions
numeroSousRegion = list()
numeroSousRegion = [x for x in range(1,nbrTotalSousRegion+1)]*int((nbrLociByRegion/nbrSousRegion))
numeroSousRegion.sort()

# recherche des primers universels souhaités
primer = [primerUniv[x] for x in primerUniv.keys() if x == PrimerU]
primer = primer[0]

# Remplissage de la classe LocusDataClass avec tous les Loci nécécesaires
total_locus = list()
for sousReg, reg, i, start, end, code in zip(numeroSousRegion,numeroRegion,range(nbrLoci),startPositions,endPositions,concatenate_CodeLocus):
  total_locus.append(LocusDataClass(sousRegionN=sousReg, regionN=reg, locusN=i+1, chrName=chromosome, startSeq=start, endSeq=end, primers_Univ = primer,codeLocus=code))

# Attribution du code combinatoire de la sous-région pour chaque locus 
for code, locus in zip(concatenate_CodeSousRegion, total_locus) :
    locus.codeSousRegion = code

# Attribution des sequences d'ADN complémentaires sondes primaires par locus en fonction des coordonnées génomiques
for locus in total_locus :
    temp = []
    for seq in listSeqGenomic :
            if locus.startSeq <= int(seq[0]) and int(seq[1])< locus.endSeq:
                temp.append([seq[0],seq[2]])                    
            else :
                pass
    random.shuffle(temp)
    temp = temp[:nbrProbeByLocus]
    temp.sort()
    locus.seqProbe = [x[1] for x in temp[:nbrProbeByLocus]]

# Affichage pour exemple d'un locus :
locus = total_locus[0].__dict__
[print(x,':',locus[x]) for x in locus.keys()]

#%%COMPLETION DES SEQUENCES PRIMAIRES AVEC CODES LOCUS,CODES REGIONS ET PRIMERS UNIVERSELS

import copy
# Insertion des binding sites pour les barcodes des loci selon le schéma suivant :
# Bcdx_Bcdy_SeqADNgenomic_Bcdx_Bcdy
for locus in total_locus :
    count = 0
    tempSeq = []
    tempBcd = []
    for i in locus.codeLocus : #000010000100000000000000000000'
        if i == '0' :
            count += 1
        elif i== '1' :
            tempBcd.append(barcodes[count][0])
            tempSeq.append(barcodes[count][1])
            count +=1
    locus.bcdLocus = tempBcd
    seqWithBcd = []
    for item in locus.seqProbe :
        seqWithBcd.append(tempSeq[0]+' '+tempSeq[1]+' '+ item +' '+tempSeq[0]+' '+tempSeq[1])
    locus.seqProbe = seqWithBcd

# Insertion de la sequence barcode pour la région selon le schéma suivant:
# Bcd-region_(Bcdx_Bcdy_SeqADNgenomic_Bcdx_Bcdy)
# ATTENTION : Il faut 2 Bcds pour identifier la région en marquage combinatoire,
# La moitié des sondes primaires seront marquées pour l'un ou pour l'autre.
bitsUseForBarcodeCombi = len(codeBookLocus[0])
for locus in total_locus :
    count = 0+bitsUseForBarcodeCombi #on utilisprint("exemple d'une séquence primaire :")e les Bcds qui ne sont pas pris pour le combinatoire des Loci
    tempSeq = []
    tempBcd = []
    for i in locus.codeSousRegion :
        if i == '0' :
            count += 1
        elif i== '1' :
            tempBcd.append(barcodes[count][0])
            tempSeq.append(barcodes[count][1])
            count +=1
    locus.bcdRegion = tempBcd
    seqWithBcdReg = []
    countSeq = 1
    for item in locus.seqProbe :
        if countSeq % 2 == 0 :
            seqWithBcdReg.append(tempSeq[0]+' '+ item)
            countSeq +=1
        else :
            seqWithBcdReg.append(tempSeq[1]+' '+ item)
            countSeq +=1
    locus.seqProbe = seqWithBcdReg

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
resultDetails = rootFolder+os.sep+'1_Library_details'
with open (resultDetails, 'w') as file :
    for locus in total_locus :
        file.write('Chromosome:'+str(locus.chrName)+' Locus_N°'+str(locus.locusN)\
+' Start:'+str(locus.startSeq)+' End:'+str(locus.endSeq)+' Bcd_locus:'+locus.bcdLocus[0]\
+'/'+locus.bcdLocus[1]+' Region_SousRegion_N°'+str(locus.regionN)+'_'+str(locus.sousRegionN)\
+' Bcd_Sous_Region:'+locus.bcdRegion[0]+'/'+locus.bcdRegion[1]+'\n')
        for seq in locus.seqProbe :
            file.write(seq+'\n')
            
#fichier avec toutes les séquences (sans espace) uniquement : Full_sequence_Only
fullSequence = rootFolder+os.sep+'2_Full_sequence_Only'
with open (fullSequence, 'w') as file :
    for locus in total_locus :
        for seq in locus.seqProbe :
            file.write(seq.replace(' ','')+'\n')
            
#fichier avec résumé des informations (sans séquence) : Library_Summary
Summary = rootFolder+os.sep+'3_Library_Summary.csv'
with open (Summary, 'w') as file :
    file.write('Chromosome,Locus_N°,Start,End,Code_locus,Bcd_locus,\
Region_SousRegion_N°,Code_Region,Bcd_Sous_Region,PU.Fw,PU.Rev\n')
    for locus in total_locus :
        file.write(str(locus.chrName)+','+str(locus.locusN)+','+str(locus.startSeq)\
+','+str(locus.endSeq)+','+"'"+str(locus.codeLocus)+','+locus.bcdLocus[0]+'/'+locus.bcdLocus[1]\
+','+str(locus.regionN)+'_'+str(locus.sousRegionN)+','+"'"+str(locus.codeSousRegion)\
+','+locus.bcdRegion[0]+'/'+locus.bcdRegion[1]+','+locus.primers_Univ[0]+','\
+locus.primers_Univ[2]+'\n')  
