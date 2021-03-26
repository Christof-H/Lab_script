#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  9 13:55:57 2021

@author: julian
"""

import numpy as np
# from functionsForLateAnalysis import bootstrap, plotMeanMatrix, plotMatrixSelectedBarcodes, calculateKDE, calculateAllKDEs
# from functionsForLateAnalysis import KS, power_law, hox_kiss, pearson_correlation, bootstrap2
import matplotlib.pyplot as plt
import os
import json
from datetime import datetime
import argparse
import csv
#from matrixOperations.alignBarcodesMasks import plotMatrix
from scipy.io import loadmat
import pandas as pd
import seaborn as sns
import scipy
import pylab
from scipy.optimize import curve_fit
import copy

path = '/mnt/PALM_dataserv/DATA/Christophe/Julian'
barcodes = [3,4,5,6,7,13,23,25,37,43,48,60,62,66,67,68,90,91,92]
data = loadmat(path + os.sep + 'HiMscMatrix.mat') #Change to use new directory
genomicCoords = loadmat('/mnt/PALM_dataserv/DATA/Christophe/Julian/coord_3R.mat')
genomicDistanceMatrix = np.zeros([len(barcodes),len(barcodes)])
genomicDistanceThreshold = 1.5e5

#genomicDistanceMatrix = matrice de la distance genomique entre RTs en paire de bases
for i in range(len(barcodes)):
     for j in range(0,i):                     
        bin1 = barcodes[i]-1        
        bin2 = barcodes[j]-1
        genomicDistanceMatrix[i][j] = abs(genomicCoords['coordlib3R'][bin1][1]-genomicCoords['coordlib3R'][bin2][1])
        
# np.where(+condition) retourne l'indexe en Y,X  de la cellule qui respecte la condition
indexes = np.where((genomicDistanceMatrix<=genomicDistanceThreshold) & (genomicDistanceMatrix>0))
matrix1 = np.zeros([len(indexes[0]),data['distanceMatrixCumulative'].shape[2]])

for i in range(len(indexes[0])):
    
    bin1 = indexes[0][i]
    bin2 = indexes[1][i]
    
    matrix1[i,:] = data['distanceMatrixCumulative'][bin1,bin2,:]

matrix2 = np.reshape(matrix1,[1,matrix1.shape[0]*matrix1.shape[1]])
vectorDistances = matrix2[~np.isnan(matrix2)]

# n: is the number of counts in each bin of the histogram
# bins: is the left hand edge of each bin
# patches is the individual patches used to create the histogram, e.g a collection of rectangles

n ,bins, patches = plt.hist(vectorDistances, range=(0,1.5), bins = 30, align = 'mid' )
plt.setp(patches[0:5], 'facecolor', 'r') 
plt.xlabel('Distance, µm')
plt.ylabel('Frequence')
plt.title(f'Genomic Distance Threshold  {str(genomicDistanceThreshold)}')

#%% Calcul de la frequence des loci à une distance inférieure ou égal a 250nm
total = 0
for i in range(len(n)):
    total = total + int(n[i])
print(total)
inf250nm = 0
for i in range(5) :
    inf250nm = inf250nm + int(n[i])
print(inf250nm)
pourcentage = round(inf250nm*100/total,0)
print (f"pourcentage : {pourcentage}%")