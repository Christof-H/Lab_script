#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 11 17:07:29 2021

@author: christophe
"""
import numpy as np
import os

locus_center=np.array([[1,10969843],[2,13034877],[3,15664949],[4,17894946]])

# Création d'un array rempli de 0 de la taille du tableau désiré
genomic_distance = np.zeros([4,4])

for i in range(0,4):
    for j in range(0,4):
        genomic_distance[i,j] = abs(locus_center[i,-1]-locus_center[j,-1])
print(genomic_distance)

rootFolder= '/home/christophe/Documents/Informatique/Python/Scripts/Genomic_distance'
file_path = rootFolder + os.sep + 'genomic_distance_array'
np.save(file_path,genomic_distance,'arr') # arr = array_like, Array data to be saved


# Pour charger l'array :
# Vérifier que vous êtes bien dans le bon dossier avec pwd et lancer :
# genomic_distance = np.load('genomic_distance.npy')