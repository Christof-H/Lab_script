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
1-Le script permet de générer tous les codes de bits à partir d'un nombre 
défini de bits par l'utilisateur.
2-Le script tri les bits en fonction d'un Hamming weight défini par l'utilisateur
3-L'utilisateur détermine le nombre de code de bits qu'il veut utiliser par rapport
à la liste totale de code de bits disponible
4-le script va recherche par force brute, un ordre de code de bits tout en calculant 
un score de bits en commun avec les codes adjacents (nbre de codes adjacents défini 
par l'utilisateur) afin d'obtenir la meilleure séquence de code avec le moins de
bits en commun 
5-Le script va essayer de remplacer les codes qui possèdent des bits en commun,
avec des codes de bits qui n'avaient pas été utlisisés pour générer la première
séquence de codes.

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

import numpy as np 
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

print ('-'*70)
list_totalbits = Code_generation(bits,Nbr_weight)

#%%------------BLOC 2 : Obtention de la meilleure séquence de code-------------

# Génère un ordre aléatoire de code de bits et calcul un Hamming score, ceci
# répété un certain nombre de fois et enregiste la séquence de code ou il y a
# le plus petit Hamming score 



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

nbrLoci = int(input("De combien de codes avez vous besoin (ex:10 loci = 10 codes) ?\n:"))
distance_inter_bit = int(input("quelle est le nombre de codes successifs qui ne doivent \
pas présenter le même bit ?\n:"))

# Nombre de fois ou l'on va tester des combinaisons de codes aléatoires
# number_trial = range(1,xxx) xxx = nbre de tests
number_trial = range(1,100)

best_bits_around = 10000000
best_list_bits_order = []

for trial in number_trial :
    random.shuffle(list_totalbits)
    partial_list_totalbits = deepcopy(list_totalbits[:nbrLoci])
    reste_list_totalbits =  deepcopy(list_totalbits[nbrLoci:])
    i1 = 0
    total_commun_bits = 0
    commun_bits = 0
    while i1 < len(partial_list_totalbits):
        if i1 < distance_inter_bit :
            nbrRange =distance_inter_bit + 1
            commun_bits = CheckBitRangeTwoWays (1, nbrRange, partial_list_totalbits[i1],partial_list_totalbits,i1)
            total_commun_bits = total_commun_bits + commun_bits
            # print(commun_bits)
            i1 +=1
        elif i1 >= distance_inter_bit and i1 < (len(partial_list_totalbits)-distance_inter_bit):
            nbrRange =distance_inter_bit + 1
            commun_bits = CheckBitRangeOneWay (1,nbrRange,partial_list_totalbits,i1)
            total_commun_bits = total_commun_bits + commun_bits
            # print(commun_bits)
            i1 +=1
        elif i1 >= (len(partial_list_totalbits)-distance_inter_bit) and i1 != (len(partial_list_totalbits)-1):
            nbrRange = len(partial_list_totalbits)-i1
            commun_bits = CheckBitRangeOneWay (1,nbrRange,partial_list_totalbits,i1)
            total_commun_bits = total_commun_bits + commun_bits
            # print(commun_bits)
            i1 +=1
        elif i1 == (len(partial_list_totalbits)-1) :
            count = 0
            # print(count)
            i1 +=1
    print (f"{trial} - Total des bits en commun sur toutes les comparaisons : {total_commun_bits}")
    if total_commun_bits < best_bits_around :
        best_bits_around = deepcopy(total_commun_bits)
        best_list_bits_order = deepcopy(partial_list_totalbits)
        reste_from_best_list_bits_order = deepcopy(reste_list_totalbits)
print('-'*70)
print ("La séquence des codes ayant obtenu le meilleur résultat est stocké \
dans la variable : best_list_bits_order")
print ("meilleur total des bits :", best_bits_around)
print (f"Avec {bits} bits utilisés, pour {nbrLoci} codes utilisés \n\
et une distance entre des bits en commun d'au moins {distance_inter_bit} codes")
#%%-------BLOC 3 : Affichage des bits en commun -5 à +5 pour les 1er bits------

# Calcul du score de bits en commun en vérifiant les premiers avec les derniers
# comme si le codes de bits était l'un derrière l'autre sans gap

# distance_inter_bit = int(input("quelle est le nombre de codes successifs qui ne doivent \
# pas présenter le même bit\n:"))
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
    print (f"Avec {bits} bits utilisés, pour {nbrLoci} codes utilisés \n\
et une distance entre des bits en commun d'au moins {distance_inter_bit} codes")

Affichage_bits_commun (best_list_bits_order)

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
print (f"Avec {bits} bits utilisés, pour {nbrLoci} codes utilisés\n\
et une distance entre des bits en commun d'au moins {distance_inter_bit} codes")


#%%------BLOC 5 : Remplacement des codes qui des bits en commun----------------

# A partir de l'odre de bit qui a obtenu le meilleur score. On remplacle les bits
# avec ceux que nous n'avons pas utilisé pour générer le premier ordre (reste_from_best_list_bits_order)
     
  
# Fonction qui check les codes de bits qui ont des bits en commun et qui
# remplace le code par un bit de la liste de bits non utilisés
def Replace_bit (listeResteBits,ordered_list, nbrRange) :
    random.shuffle(listeResteBits)
    y=0
    while y < len(listeResteBits):
        count_newbit = CheckBitRangeTwoWays (1,nbrRange,listeResteBits[y],ordered_list,i)
        if count_newbit != 0 :
            y +=1
            find_bit = False
        elif count_newbit == 0 :
            print(f"Replacement de {ordered_list[i]} par {listeResteBits[y]}",'i=',i)
            ordered_list[i] ,listeResteBits[y] = listeResteBits[y], ordered_list[i]
            find_bit = True
            break
    if find_bit == False :
        print('Auncun code pour remplacer', 'i=',i)
    else :
        pass       


# distance_inter_bit = int(input("quelle est le nombre de codes successifs qui ne doivent \
# pas présenter le même bit\n:"))

            
i=0            
while i < len(best_list_bits_order):
     # Condition pour les codes : best_list_bits_order[0] à best_list_bits_order[distance_inter_bit]        
    if i < distance_inter_bit :
        nbrRange =distance_inter_bit + 1
        commun_bits = CheckBitRangeTwoWays (1, nbrRange, best_list_bits_order[i],best_list_bits_order,i)
        if commun_bits == 0:
            print('0 bit en commun, pas de remplacement', 'i=', i)
            i +=1
        if commun_bits != 0 :
            Replace_bit (reste_from_best_list_bits_order,best_list_bits_order, nbrRange)
            i +=1
# Condition pour les codes jusqu'à best_list_bits_order[liste totale - distance_inter_bit]        
    elif i >= distance_inter_bit and i < (len(best_list_bits_order)-distance_inter_bit):
        nbrRange =distance_inter_bit + 1
        commun_bits = CheckBitRangeOneWay (1, nbrRange, best_list_bits_order,i)
        if commun_bits == 0:
            print('0 bit en commun,pas de remplacement', 'i=', i)
            i +=1
        if commun_bits != 0 :
            Replace_bit(reste_from_best_list_bits_order,best_list_bits_order,nbrRange)
            i +=1
# condition pour les bits de fin
    elif i >= (len(best_list_bits_order)-distance_inter_bit) and i != (len(best_list_bits_order)-1):
        nbrRange = len(best_list_bits_order)-i
        commun_bits = CheckBitRangeOneWay (1,nbrRange,best_list_bits_order,i)
        if commun_bits == 0:
            print('0 bit en commun,pas de remplacement', 'i=', i)
            i +=1
        if commun_bits != 0 :
            Replace_bit(reste_from_best_list_bits_order,best_list_bits_order,nbrRange)
            i +=1
    elif i == (len(best_list_bits_order)-1) :
        print('0 bit en commun,pas de remplacement','i=',i)
        i +=1
print('-'*70)
print(f"Avertissement : refaire bloc 3 si des codes ont été changé dans les \
{distance_inter_bit} premiers et derniers codes.")



rep = input("Voulez-vous relancer :\n\
1-un affichage des bits en commun pour la nouvelle séquence de bits (taper : 1)\n\
2-Sortir (taper : 2)\n\
Réponse (1 ou 2) :")
if rep == '1':
    Affichage_bits_commun (best_list_bits_order)    
else :
    pass




