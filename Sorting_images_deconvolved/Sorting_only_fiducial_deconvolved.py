#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  6 15:18:53 2022

@author: christophe

This script is based solely on the fiducial images of the different cycles.
It will allow a segmentation of the fiducial spots as for a barcode. 
We will thus be able to check if we detect the same number of fiducial spots 
in the first cycles as in the last from the segmentedObjects file.
The script will :
- create links of the fiducial images (scan_RTX_00X_ROI_converyed_decon_ch00.tif) 
in folders ROI001, ROI002...
- duplicate these fiducial images and rename them, to identify them as barcodes (scan_...._converyed_decon_ch00.tif 
to scan_...._converyed_decon_ch01.tif)



"""

import os, glob

rootFolderRTs = '/mnt/grey/DATA/ProcessedData_2022/Experiment_13_Tof_tests/RTs_deconv_other_exp'
rootFolderDapi = '/mnt/grey/DATA/ProcessedData_2022/Experiment_13_Tof_tests/Dapi_deconv_other_exp'
rootFolder='/mnt/grey/DATA/AnalyzedData_2022/Experiment_13_Tof_tests'

destFolder = rootFolder + os.sep + 'Images_deconvolved_sorted'
if not os.path.isdir(destFolder):
    os.mkdir(destFolder)
    
#%%TRI TOUS LES RTs EN FONCTION DES TOUS LES ROIS 
#Trie des images en fonctions du ROI : classe les images de tous les RTs d'un
# même ROI dans un répertoire
destFolderfiducial_ch00_and_ch01 = destFolder+ os.sep +'fiducial_ch00_and_ch01'
if not os.path.isdir(destFolderfiducial_ch00_and_ch01):
    os.mkdir(destFolderfiducial_ch00_and_ch01)
# Regroupement des noms des images DAPI et RTs dans une liste filesDAPI et filesRT
# filesDAPI = ['mnt/.../Dapi_deconvolved/scan_007_DAPI_022_ROI_converted_decon_ch01.tif',
# '/mnt/.../Dapi_deconvolved/scan_007_DAPI_025_ROI_converted_decon_ch00.tif',.......]
#
# filesRT = ['/mnt/.../RT_deconvolved/scan_008_RT_3_4_022_ROI_converted_decon_ch00.tif',
# '/mnt/.../RT_deconvolved/scan_008_RT_11_12_025_ROI_converted_decon_ch01.tif',.......]
filesDAPI = glob.glob(rootFolderDapi + os.sep + "*.tif")
filesFiducial = glob.glob(rootFolderRTs + os.sep + "*ch00.tif")

# Création d'une liste de tous les numéro de ROI en faisant un split sur les noms contenus 
# dans filesDAPI
# listROIs = ['022', '025', '014'....]
listeROIs = []
for x in filesDAPI:
    ROI = os.path.basename(x).split('_')[6]
    if ROI not in listeROIs:
        listeROIs.append(ROI)

# Création des dossiers correspondants aux différents ROI : ROI001, ROI002:
# os.path.isdir(path): Return True if path is an existing directory.
for ROI in listeROIs :
    newFolder = destFolderfiducial_ch00_and_ch01 + os.sep + "ROI" + ROI
    if not os.path.isdir(newFolder):
        os.mkdir(newFolder)
# Trie des images Dapi
    filesDAPIFiltered = [x for x in filesDAPI if os.path.basename(x).split('_')[3]== ROI]
    for x in filesDAPIFiltered:
        y = newFolder + os.sep + os.path.basename(x)
        os.system('ln -s ' + x + " " + y)
# Trie des images RTs
    filesFiltered = [x for x in filesFiducial if os.path.basename(x).split('_')[3] == ROI]
    for x in filesFiltered:
        y = newFolder + os.sep + os.path.basename(x)
        z = newFolder + os.sep + os.path.basename(x).replace("ch00", "ch01")
        os.system('ln -s ' + x + " " + y)
        os.system('ln -s ' + x + " " + z)