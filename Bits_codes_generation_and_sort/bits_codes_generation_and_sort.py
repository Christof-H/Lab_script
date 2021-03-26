#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 19 16:27:52 2021

@author: christophe

Ce script génère un code de bits avec le moins de bits en commun avec les codes de bits adjacents.
--------------------------------------------------------------------
ATTENTION : ce scipt recherche la meilleure combinaison de codes comme si cette séquence
de codes de bits était collée l'une à la suite de l'autre.
Les premiers codes vont être comparés avec les derniers 
--------------------------------------------------------------------
BLOC 1-Le script permet de générer tous les codes de bits à partir d'un nombre 
défini de bits par l'utilisateur, et ne retient les codes ne présentant qu'un certains
nombre de bit (défini par le Hamming weight)

BLOC 2-Le script choisi au hasard un nombre de codes (défini par l'utilisateur),
va tester la meilleure combinaison et essayer de remplacer les codes qui possèdent des
bits en commun avec les codes adjacents. Si il n'y arrive pas, le script génère 
un nouveau groupe de codes au hasard et refait de même. 

Au final vous obtenez la meilleure séquence de codes dans la variable : best_list_bits_order

"""
# Python Bitwise Operators
# Avant et après les opétateurs, il faudra placer des chiffres en binaire
# example 1=1, 2=10, 3=11, 4=100, 5=101..... chiffre en base binaire = nombre
# le chiffre binaire de 4 est '100'. bin(4) = 100
# A l'inverse pour trouver le chiffre de base binaire à partir du nombre
# on utilise la fonction int('string',base) base binaire = 2
# int('100',2) = 4
# & : Operator copies a bit to the result if it exists in both operands
# example : si on veut comparer 1011 (binaire 11) à 1101 (binaire 13):
# 11 & 13 = 9 bin(9) = 1001 ou int('1011',2)&int('1101',2) = 9
# ^ : It copies the bit if it is set in one operand but not both.
# example : 11^13=6, bin(6) = 110 ou int('1011',2)^int('1101',2) = 6




# ------BLOC 1 : Génère et trie les codes en fonction du Hamming Weight--------

import os
import random
from copy import deepcopy

# Fonction qui permet de comparer 2 bit codes pour révéler les bits en commun
def NbrBitCommun (a,b):
    result = 0
    similar_bits = int(a,2) & int(b,2) #recherche et garde les bits en commun ex 1001 & 1100 = 1000 soit 8 car bin(8) = 0b1000
    for i in bin(similar_bits)[2:]:    # on compte le nombre de 1 après 0b ex : 0b1000 soit result = 1
        result = result + int(i)
    return int(result)

# Fonction qui permet de générer tous les codes à partir d'un nombre de bits, et
# qui trie les codes en fonction d'un Hamming weight
def Code_generation(n, hamming_weight):
    result = []
    total = []
    count_total = 0
    for i in range(1<<n):
        s=bin(i)[2:]
        s='0'*(n-len(s))+s
        # total.append(s)  #activer cette ligne si vous voulez garder tous les bits possible ATTENTION mémoire
        count_total += 1
        count = 0
        for i in s :
            if i == '1' :
                count += 1
        if count == hamming_weight :
            result.append(s)
    print(f"Pour {n} bits, il y a au total {count_total} codes générés.\n\
Après un tri avec un hamming weight de {hamming_weight}, il reste {len(result)} codes.")
    return result

# Paramétrage du nombre de bits et du Hamming codes désiré
bits = int(input("Avec combien de bits, voulez-vous générer vos codes ?\n:"))
Nbr_weight = int(input("quel est le Hamming weight que vous désirez ?\n: "))

###
print ('-'*70)
list_totalbits = Code_generation(bits,Nbr_weight)
print ('-'*70)
###


#%%-------BLOC 2 : Test et remplacement des codes présentant des bits en commun
# avec les codes adjacents

# Fonction qui calcul un score de bit en commun pour les bits situé de la
# position -5 à +5 par defaut:
def CheckBitRangeOneWay (number1,number2,ordered_list,i):
    count = 0
    for x in range(number1,number2) :
        count = count + NbrBitCommun (ordered_list[i],ordered_list[i+x])
    return count

# Fonction qui calcul un score de bit en commun pour le bit de remplacement 
# avec les bits situé de la position -x à +x :
def CheckBitRangeTwoWays (number1,number2,newbit,ordered_list,i):
    count = 0
    for x in range(number1,number2) :
        count = count + NbrBitCommun (newbit,ordered_list[i+x])
        count = count + NbrBitCommun (newbit,ordered_list[i-x])
    return count

number_Loci = int(input("De combien de codes avez vous besoin (ex:10 loci = 10 codes) ?\n:"))
distance_inter_bit = int(input("quelle est le nombre de codes successifs qui ne doivent \
pas présenter le même bit ?\n:"))

# Nombre de fois ou l'on va tester des combinaisons de codes aléatoires
# number_trial = range(1,xxx) xxx = nbre de tests
number_trial = range(1,300)

def Search_and_replace (nbr_trial, nbrLoci, dist_inter_bit, list_bits):
    for trial in nbr_trial :
        best_bits_around = 10000000
        random.shuffle(list_bits)
        partial_list_bits = deepcopy(list_bits[:nbrLoci])
        reste_list_bits =  deepcopy(list_bits[nbrLoci:])
        i1 = 0
        total_commun_bits = 0
        commun_bits = 0
        while i1 < len(partial_list_bits):
            if i1 < dist_inter_bit :
                nbrRange =dist_inter_bit + 1
                commun_bits = CheckBitRangeTwoWays (1, nbrRange, partial_list_bits[i1],partial_list_bits,i1)
                total_commun_bits = total_commun_bits + commun_bits
                i1 +=1
            elif i1 >= dist_inter_bit and i1 < (len(partial_list_bits)-dist_inter_bit):
                nbrRange =dist_inter_bit + 1
                commun_bits = CheckBitRangeOneWay (1,nbrRange,partial_list_bits,i1)
                total_commun_bits = total_commun_bits + commun_bits
                i1 +=1
            elif i1 >= (len(partial_list_bits)-dist_inter_bit) and i1 != (len(partial_list_bits)-1):
                nbrRange = len(partial_list_bits)-i1
                commun_bits = CheckBitRangeOneWay (1,nbrRange,partial_list_bits,i1)
                total_commun_bits = total_commun_bits + commun_bits
                i1 +=1
            elif i1 == (len(partial_list_bits)-1) :
                count = 0
                i1 +=1
        if total_commun_bits < best_bits_around :
            best_bits_around = deepcopy(total_commun_bits)
            best_list_bits = deepcopy(partial_list_bits)
            reste_from_best_list_bits = deepcopy(reste_list_bits)
    return (best_list_bits, reste_from_best_list_bits)





def Count_bits_commun (list_bits):
    item = 0
    total_commun_bits = 0
    commun_bits = 0
    while item < len(list_bits):
        if item < distance_inter_bit :
            nbrRange =distance_inter_bit + 1
            commun_bits = CheckBitRangeTwoWays (1, nbrRange, list_bits[item],list_bits,item)
            total_commun_bits = total_commun_bits + commun_bits
            item +=1
        elif item >= distance_inter_bit and item < (len(list_bits)-distance_inter_bit):
            nbrRange =distance_inter_bit + 1
            commun_bits = CheckBitRangeOneWay (1,nbrRange,list_bits,item)
            total_commun_bits = total_commun_bits + commun_bits
            item +=1
        elif item >= (len(list_bits)-distance_inter_bit) and item != (len(list_bits)-1):
            nbrRange = len(list_bits)-item
            commun_bits = CheckBitRangeOneWay (1,nbrRange,list_bits,item)
            total_commun_bits = total_commun_bits + commun_bits
            item +=1
        elif item == (len(list_bits)-1) :
            count = 0
            item +=1
    return (total_commun_bits)





# Calcul du score de bits en commun en vérifiant les premiers avec les derniers
# comme si le codes de bits était l'un derrière l'autre sans gap
def Affichage_bits_commun (list_bits):
    item = 0
    total_commun_bits = 0
    commun_bits = 0
    while item < len(list_bits):
        if item < distance_inter_bit :
            nbrRange =distance_inter_bit + 1
            commun_bits = CheckBitRangeTwoWays (1, nbrRange, list_bits[item],list_bits,item)
            total_commun_bits = total_commun_bits + commun_bits
            print(commun_bits)
            item +=1
        elif item >= distance_inter_bit and item < (len(list_bits)-distance_inter_bit):
            nbrRange =distance_inter_bit + 1
            commun_bits = CheckBitRangeOneWay (1,nbrRange,list_bits,item)
            total_commun_bits = total_commun_bits + commun_bits
            print(commun_bits)
            item +=1
        elif item >= (len(list_bits)-distance_inter_bit) and item != (len(list_bits)-1):
            nbrRange = len(list_bits)-item
            commun_bits = CheckBitRangeOneWay (1,nbrRange,list_bits,item)
            total_commun_bits = total_commun_bits + commun_bits
            print(commun_bits)
            item +=1
        elif item == (len(list_bits)-1) :
            count = 0
            print(count)
            item +=1
    print ('-'*70)
    print ("Score du nombre de bits trop proche pour une séquence de code : ", total_commun_bits)
    print (f"Avec {bits} bits utilisés, pour {number_Loci} codes utilisés \n\
et une distance entre des bits en commun d'au moins {distance_inter_bit} codes")
    print ("La meilleure séquence de code de bits est stockée dans la variable : best_list_bits_order")
    return (total_commun_bits)




# A partir de l'ordre de bit qui a obtenu le meilleur score. On remplace les bits
# avec ceux que nous n'avons pas utilisé pour générer le premier ordre (reste_from_best_list_bits_order)
     

# Fonction qui check les codes de bits qui ont des bits en commun et qui
# remplace le code par un bit de la liste de bits non utilisés
def Replace_bit (listeResteBits,ordered_list, nbrRange,i) :
    random.shuffle(listeResteBits)
    y=0
    while y < len(listeResteBits):
        count_newbit = CheckBitRangeTwoWays (1,nbrRange,listeResteBits[y],ordered_list,i)
        if count_newbit != 0 :
            y +=1
            find_bit = False
        elif count_newbit == 0 :
            # print(f"Replacement de {ordered_list[i]} par {listeResteBits[y]}",'i=',i)
            ordered_list[i] ,listeResteBits[y] = listeResteBits[y], ordered_list[i]
            find_bit = True
            break
    if find_bit == False :
        # print('Auncun code pour remplacer', 'i=',i)
        pass
    else :
        pass       


def Test_and_replace (list_bits_order) :           
    i=0            
    while i < len(list_bits_order):
         # Condition pour les codes : list_bits_order[0] à list_bits_order[distance_inter_bit]        
        if i < distance_inter_bit :
            nbrRange =distance_inter_bit + 1
            count_bits = CheckBitRangeTwoWays (1, nbrRange, list_bits_order[i],list_bits_order,i)
            if count_bits == 0:
                # print('0 bit en commun, pas de remplacement', 'i=', i)
                i +=1
            if count_bits != 0 :
                Replace_bit (reste_from_best_list_bits_order,list_bits_order, nbrRange,i)
                i +=1
    # Condition pour les codes jusqu'à list_bits_order[liste totale - distance_inter_bit]        
        elif i >= distance_inter_bit and i < (len(list_bits_order)-distance_inter_bit):
            nbrRange =distance_inter_bit + 1
            count_bits = CheckBitRangeOneWay (1, nbrRange, list_bits_order,i)
            if count_bits == 0:
                # print('0 bit en commun,pas de remplacement', 'i=', i)
                i +=1
            if count_bits != 0 :
                Replace_bit(reste_from_best_list_bits_order,list_bits_order,nbrRange,i)
                i +=1
    # condition pour les bits de fin
        elif i >= (len(list_bits_order)-distance_inter_bit) and i != (len(list_bits_order)-1):
            nbrRange = len(list_bits_order)-i
            count_bits = CheckBitRangeOneWay (1,nbrRange,list_bits_order,i)
            if count_bits == 0:
                # print('0 bit en commun,pas de remplacement', 'i=', i)
                i +=1
            if count_bits != 0 :
                Replace_bit(reste_from_best_list_bits_order,list_bits_order,nbrRange,i)
                i +=1
        elif i == (len(list_bits_order)-1) :
            # print('0 bit en commun,pas de remplacement','i=',i)
            i +=1
#     print('-'*70)
#     print(f"Avertissement : refaire bloc 3 si des codes ont été changé dans les \
# {distance_inter_bit} premiers et derniers codes.")

count_trial = 0
result = False
for trail in range(1,100) :
    if result == False :
        best_list_bits_order, reste_from_best_list_bits_order = Search_and_replace (number_trial, number_Loci, distance_inter_bit, list_totalbits)
        count_trial +=1
        print(f"Trial n°{count_trial}")
        for trial2 in range(1,6):
            count_bits = Count_bits_commun (best_list_bits_order)
            if count_bits !=0:
                Test_and_replace (best_list_bits_order)
            elif count_bits == 0:
                Affichage_bits_commun (best_list_bits_order)
                result = True
                break
        if result != True :
            commun_bits_test=Count_bits_commun (best_list_bits_order)
            print(f"Impossible to replace all bits....Best score for this trial = {commun_bits_test}")
            pass
    elif result != False :
        pass
#%% Enregistrement de la meilleure séquence de codes :

# Renseigner le chemin d'accès du dossier dans lequel vous voulez sauvegarder
# votre fichier contenant tous les codes de bits :        
root_folder = '/home/christophe/Documents/Informatique/Python/Scripts/Bits_code_generation_Resultats'

# Génère un nom de fichier en fonction des paramètres 
file_name = 'Best_Seq_Codes_'+str(bits)+'Bit_'+str(Nbr_weight)+'HWeiht_'+str(number_Loci)+'Cod_Esp'+str(distance_inter_bit)

fullFileName=root_folder + os.sep + file_name


# Vérification de l'existence du répertoire pour sauvegarder le fichier  
check_path = os.path.exists(root_folder)
if check_path != True :
    print ('-'*70)
    print ("The path to save the file is not correct.")
    print ("Please indicate the correct path in the variable 'root_folder'.")
    print ('-'*70)
elif check_path == True:
    fichier = open(fullFileName, "w")
    fichier.write('Séquence de '+str(number_Loci)+' codes générés avec '+str(bits)+
    ' bits, pour un Hamming weight de '+str(Nbr_weight)+' , avec un espacement de '
    +str(distance_inter_bit)+' codes entre chaque bit en commun\n')
    fichier.write('best_list_bits_order=\n')
    for item in best_list_bits_order :
        fichier.write('\''+item+'\''+','+'\n')
    fichier.close()

#%%-------------BLOC 4 : Affichage des bits en commun +1 à +5------------------ 

# Partie pour tester et affiché la combinaison de codes la plus optimale 
# trouvée précédemment sans comparer le début et la fin de la séquence de code de bits

# distance_inter_bit = int(input("quelle est le nombre de codes successifs qui ne doivent \
# pas présenter le même bit\n:"))

i2 = 0
total_commun_bits = 0
commun_bits = 0
while i2 < len(best_list_bits_order):
    if i2 < (len(best_list_bits_order)-distance_inter_bit) :
        nbrRange =distance_inter_bit + 1
        commun_bits = CheckBitRangeOneWay (1, nbrRange, best_list_bits_order,i2)
        total_commun_bits = total_commun_bits + commun_bits
        print(commun_bits)
        i2 +=1
    elif i2 >= (len(best_list_bits_order)-distance_inter_bit) and i2 != (len(best_list_bits_order)-1):
        nbrRange = len(best_list_bits_order)-i2
        commun_bits = CheckBitRangeOneWay (1,nbrRange,best_list_bits_order,i2)
        total_commun_bits = total_commun_bits + commun_bits
        print(commun_bits)
        i2 +=1
    elif i2 == (len(best_list_bits_order)-1) :
        count = 0
        print(count)
        i2 +=1
print ('-'*70)
print ("Score du nombre de bits trop proche pour une séquence de code : ", total_commun_bits)
print (f"Avec {bits} bits utilisés, pour {number_Loci} codes utilisés\n\
et une distance entre des bits en commun d'au moins {distance_inter_bit} codes")



