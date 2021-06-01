#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 26 10:53:50 2021

@author: christophe

Script qui permet de vérifier la longueur des séquences dans un fichier 
contenant toutes les séquences d'une biblliothèque (ex : 2_Full_sequence_Only')
notamment lors d'une compilation de bibliothèques différentes'

"""
import os




rootFolder = '/home/christophe/Documents/Informatique/Python/Scripts/En cours/Verif_length_and_adjust'
fullSeqOnlyFile = '2_Full_sequence_Only'
fullSeqOnly = rootFolder + os.sep + fullSeqOnlyFile

fullSeqOnlyFileComp = '2_Full_sequence_Only_Compl'
fullSeqOnlyComp = rootFolder + os.sep + fullSeqOnlyFileComp




listFullSeqOnly = []
with open(fullSeqOnly, 'r') as file:
    for line in file :
        listFullSeqOnly.append(line.replace('\n',''))


# Fonction permettant l'évaluation de la longueur des sondes primaires de toute la librairie 
def Check_Length_Seq_Diff(librairy):
    minimal_length = 500
    maximal_length = 100
    for seq in librairy:
        if len(seq) < minimal_length :
            minimal_length = len(seq)
        elif len(seq) > maximal_length :
            maximal_length = len(seq)
    difference_pourcentage = 100-(minimal_length*100/maximal_length)
    difference_nbre = maximal_length-minimal_length
    print(f'Length of smaller primary sequence = {str(minimal_length)}-mer')
    print(f'Length of bigger primary sequence = {str(maximal_length)}-mer')
    print(f'Difference in percentage : {int(difference_pourcentage)}%')
    print(f'Difference in nucleotides : {difference_nbre}-mer')
    return int(difference_pourcentage), maximal_length

diff_Pourc, max_Length = Check_Length_Seq_Diff(listFullSeqOnly)


#%%
# Fonction permettant le complétion de nucléotides aléatoires pour des sequences dont la difference de taille est supérieure a 10%
def Completion(difference_pourcentage, Max_length, librairy, Max_diff_pourcent=10) :
# Max_diff_pourcent = pourcentage de difference maximal toléré, calculé à partir de la fonction Check_Length_Seq_Diff()
# Max_length = provient normalement de la fonction Check_Length_Seq_Diff()
# librairy = liste des sequences
    import random
    seq_completion = list()
    if difference_pourcentage >= Max_diff_pourcent :
        for seq in librairy :
            diff_seq_with_max = Max_length - len(seq)
            seq_added = ''
            for i in range (diff_seq_with_max):
                seq_added = seq_added+random.choice('atgc')
            seq_completion.append(seq+seq_added)
        print('-'*70)
        print('Completion finished')
        print('-'*70)
    else :
        print('-'*70)
        print('No completion required')
        print('-'*70)
    return seq_completion


finalFullSeqOnly = list()
finalFullSeqOnly =Completion(diff_Pourc, max_Length, listFullSeqOnly)         

#%%
# Ecriture des sequences complétées dans un nouveau fichier : 2_Full_sequence_Only_Compl

with open (fullSeqOnlyComp, 'w') as file :
    for seq in finalFullSeqOnly :
        file.write(seq+'\n')
        
        