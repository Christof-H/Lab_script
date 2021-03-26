#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 15 14:35:42 2021

@author: christophe
"""

import os, glob

rootFolderRTs = '/mnt/grey/DATA/rawData_2021/Experiment_13_Tof/RT/RT_Deconvolved'
rootFolderDapi = '/mnt/grey/DATA/rawData_2021/Experiment_13_Tof/DAPI/Dapi_deconvolved'
rootFolder='/mnt/grey/DATA/rawData_2021/Experiment_13_Tof'

destFolder = rootFolder + os.sep + 'Images_deconvolved_sorted'
os.mkdir(destFolder)

#%% Trie des images en fonctions du ROI : classe les images de tous les RTs d'un
# même ROI dans un répertoire
destFolderSortedByROI = destFolder+ os.sep +'sorted_by_ROI'
os.mkdir(destFolderSortedByROI)
# Regroupement des noms des images DAPI et RTs dans une liste filesDAPI et filesRT
# filesDAPI = ['mnt/.../Dapi_deconvolved/scan_007_DAPI_022_ROI_converted_decon_ch01.tif',
# '/mnt/.../Dapi_deconvolved/scan_007_DAPI_025_ROI_converted_decon_ch00.tif',.......]
#
# filesRT = ['/mnt/.../RT_deconvolved/scan_008_RT_3_4_022_ROI_converted_decon_ch00.tif',
# '/mnt/.../RT_deconvolved/scan_008_RT_11_12_025_ROI_converted_decon_ch01.tif',.......]
filesDAPI = glob.glob(rootFolderDapi+ os.sep+ "*.tif")
filesRT = glob.glob(rootFolderRTs + os.sep+ "*.tif")

# Création d'une liste de tous les numéro de ROI en faisant un split sur les noms contenus 
# dans filesDAPI
# listROIs = ['022', '025', '014'....]
listeROIs = []
for x in filesDAPI:
    ROI = os.path.basename(x).split('_')[3]
    if ROI not in listeROIs:
        listeROIs.append(ROI)

# Création des dossiers correspondants aux différents ROI : ROI001, ROI002:
# os.path.isdir(path): Return True if path is an existing directory.
for ROI in listeROIs :
    newFolder = destFolderSortedByROI + os.sep + "ROI" + ROI
    if not os.path.isdir(newFolder):
        os.mkdir(newFolder)
# Trie des images Dapi
    filesDAPIFiltered = [x for x in filesDAPI if os.path.basename(x).split('_')[3]== ROI]
    for x in filesDAPIFiltered:
        y = newFolder + os.sep + os.path.basename(x)
        os.system('ln -s ' + x + " " + y)
# Trie des images RTs
    filesRTFiltered = [x for x in filesRT if os.path.basename(x).split('_')[3] == ROI]
    for x in filesRTFiltered:
        y = newFolder + os.sep + os.path.basename(x)
        os.system('ln -s ' + x + " " + y)

#%% Trie des images en fonction des RTs pour tous les ROIs : classe les images 
# de tous les ROIs seulement pour certains RTs dans un même dossier
destFolderSortedForRTs = destFolder+ os.sep +'All_ROI_Sorted_For_Some_RTs'
os.mkdir(destFolderSortedForRTs)
# Regroupement des noms des images DAPI et RTs dans une liste filesDAPI et filesRT
# filesDAPI = ['mnt/.../Dapi_deconvolved/scan_007_DAPI_022_ROI_converted_decon_ch01.tif',
# '/mnt/.../Dapi_deconvolved/scan_007_DAPI_025_ROI_converted_decon_ch00.tif',.......]
#
# filesRT = ['/mnt/.../RT_deconvolved/scan_008_RT_3_4_022_ROI_converted_decon_ch00.tif',
# '/mnt/.../RT_deconvolved/scan_008_RT_11_12_025_ROI_converted_decon_ch01.tif',.......]
filesDAPI = glob.glob(rootFolderDapi+ os.sep+ "*.tif")
filesRT = glob.glob(rootFolderRTs + os.sep+ "*.tif")

# Liste des RTs que l'on veut pour tous les ROIs : 
listeRTs = ['RT1', 'RT2', 'RT3', 'RT4']

for x in filesDAPI:
    y = destFolderSortedForRTs + os.sep + os.path.basename(x)
    os.system('ln -s ' + x + " " + y)
for RT in listeRTs :
    filesRTFiltered = [x for x in filesRT if os.path.basename(x).split('_')[2] == RT]
    for x in filesRTFiltered:
        y = destFolderSortedForRTs + os.sep + os.path.basename(x)
        os.system('ln -s ' + x + " " + y)    

#%% Trie des images en fonctions de certains RTs et de certains ROIs choisis :
destFolderSortedForRTsForROIs = destFolder+ os.sep +'Some_RTs_Sorted_For_Some_ROIs'
os.mkdir(destFolderSortedForRTsForROIs)
# Regroupement des noms des images DAPI et RTs dans une liste filesDAPI et filesRT
# filesDAPI = ['mnt/.../Dapi_deconvolved/scan_007_DAPI_022_ROI_converted_decon_ch01.tif',
# '/mnt/.../Dapi_deconvolved/scan_007_DAPI_025_ROI_converted_decon_ch00.tif',.......]
#
# filesRT = ['/mnt/.../RT_deconvolved/scan_008_RT_3_4_022_ROI_converted_decon_ch00.tif',
# '/mnt/.../RT_deconvolved/scan_008_RT_11_12_025_ROI_converted_decon_ch01.tif',.......]
filesDAPI = glob.glob(rootFolderDapi+ os.sep+ "*.tif")
filesRT = glob.glob(rootFolderRTs + os.sep+ "*.tif")

# Liste des RTs que l'on veut pour certains ROIs : 
listeRTs = ['RT1', 'RT2', 'RT3', 'RT4']

# Liste des ROIs que l'on veut utiliser : 
listeROIs = ['001', '002', '003', '004', '005', '010']

for ROI in listeROIs :
# Trie des images Dapi
    filesDAPIFiltered = [x for x in filesDAPI if os.path.basename(x).split('_')[3]== ROI]
    for x in filesDAPIFiltered:
        y = destFolderSortedForRTsForROIs + os.sep + os.path.basename(x)
        os.system('ln -s ' + x + " " + y)
# Trie des images RTs en fonction des ROIs
    filesRTFiltered = [x for x in filesRT if os.path.basename(x).split('_')[3] == ROI]
  # Trie des images RTs en fonction des RTs choisis  
    for RT in listeRTs :
        finalfilesRTFiltered = [x for x in filesRTFiltered if os.path.basename(x).split('_')[2] == RT]
        for x in finalfilesRTFiltered:
            y = destFolderSortedForRTsForROIs + os.sep + os.path.basename(x)
            os.system('ln -s ' + x + " " + y)
 
