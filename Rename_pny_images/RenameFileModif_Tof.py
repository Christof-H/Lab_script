#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 13:13:10 2020

@author: christophe
"""


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 28 10:55:13 2020
@author: marcnol
"""

import os, glob

def makeLink(newFolder, x, list1, ext):
    y = newFolder + os.sep + "_".join(list1)+"."+ext
    os.system('ln -s ' + x + " " + y)
    return 1


rootFolderDAPI = "/mnt/grey/DATA/rawData_2020/Exp_Combinatory_3_Tof/DAPI/Dapi_deconvolved"
rootFolderRT = "/mnt/grey/DATA/rawData_2020/Exp_Combinatory_3_Tof/RT/RT_deconvolved"
destFolder = "/mnt/grey/DATA/rawData_2020/Exp_Combinatory_3_Tof/ROI_Tri_Tof"

# Regroupement des noms des images DAPI et RTs dans une liste filesDAPI et filesRT
# filesDAPI = ['mnt/.../Dapi_deconvolved/scan_007_DAPI_022_ROI_converted_decon_ch01.tif',
# '/mnt/.../Dapi_deconvolved/scan_007_DAPI_025_ROI_converted_decon_ch00.tif',.......]
#
# filesRT = ['/mnt/.../RT_deconvolved/scan_008_RT_3_4_022_ROI_converted_decon_ch00.tif',
# '/mnt/.../RT_deconvolved/scan_008_RT_11_12_025_ROI_converted_decon_ch01.tif',.......]


filesDAPI = glob.glob(rootFolderDAPI+ os.sep+ "*.tif")
filesRT = glob.glob(rootFolderRT + os.sep+ "*.tif")

# Création d'une liste de tous les numéro de ROI en faisant un split sur les noms contenus 
# dans filesDAPI
# listROIs = ['022', '025', '014'....]
#os.path.basename()=return the base name (last name) 
listROIs = []
for x in filesDAPI:
    ROI = os.path.basename(x).split('.')[0].split('_')[3]
    if ROI not in listROIs:
        listROIs.append(ROI)
        
# Création des dossiers pour chaque ROI à partir du split du nom dans filesDAPI : 
# newfolder = /mnt/.../ROI_Tri_Tof/ROI001
# newfolder = /mnt/.../ROI_Tri_Tof/ROI002
        
for ROI in listROIs:
    
    newFolder = destFolder + os.sep + "ROI" + ROI
    if not os.path.isdir(newFolder):
        os.mkdir(newFolder)

# Indentation : pour chaque ROI in listROIs, filesDAPIFiltered = toutes les images correspondantes
# au ROI et filesRTFiltered = toutes les images correspondantes au ROI
# pour ROI 001 filesDAPIFiltered = ['/mnt/.../Dapi_deconvolved/DAPI_001_ROI_converted_decon_ch00.tif,...]
#pour ROI 001 filesRTFiltered = ['/mnt/.../RT_deconvolved/RT_9_10_001_ROI_converted_decon_ch01.tif,...]

    filesDAPIFiltered = [x for x in filesDAPI if os.path.basename(x).split('.')[0].split('_')[3] == ROI] 
    
#???    filesRTFiltered = [x for x in filesRT if os.path.basename(x).split('.')[0].split('_')[3] == ROI] 
    
    filesRTFiltered = [x for x in filesRT if os.path.basename(x).split('.')[0].split('_')[5] == ROI] 
    
    # print([os.path.basename(x) for x in filesDAPIFiltered])
    # print([os.path.basename(x) for x in filesRTFiltered])
# Indentation : pour chaque x dans filesDAPIFiltered un soft-link est créer:
# y = /mnt/.../ROI001/DAPI_001_ROI_converted_decon_ch00.tif
# création du lien : $ ln -s mnt/.../Dapi_deconvolved/DAPI_001_ROI_converted_decon_ch00.tif /mnt/.../ROI001/DAPI_001_ROI_converted_decon_ch00.tif
    
    for x in filesDAPIFiltered:
        y = newFolder + os.sep + os.path.basename(x)
        os.system('ln -s ' + x + " " + y)

# for x in filesRTFiltered (/mnt/.../RT_deconvolved/RT_9_10_001_ROI_converted_decon_ch01.tif)
# filename = 'scan_008_RT_9_10_001_ROI_converted_decon_ch01.tif'
# list0 = ['scan', '008', 'RT', '9', '10', '001', 'ROI', 'converted', 'decon', 'ch01']
# ext = 'tif'
       
    counter=0
    for x in filesRTFiltered:
        filename=os.path.basename(x)
        list0 = filename.split(".")[0].split("_")
        ext=filename.split(".")[1]
#---------------------------------------------------------------------------
# ATTENTION POUR CETTE MANIP
# scan_008_RT_1_2_001_ROI_converted_decon_ch01.tif
# scan_008_RT_1_2_001_ROI_converted_decon_ch02.tif 
# ch01=RhoRed (=RT2) ch02=Al647 (=RT1)
#---------------------------------------------------------------------------

# ex pour filename = 'scan_008_RT_9_10_001_ROI_converted_decon_ch01.tif'
# list0 = ['scan', '008', 'RT', '9', '10', '001', 'ROI', 'converted', 'decon', 'ch01']
# if list0[-1]=="ch01"
# list1 = ['scan', '008', 'RT10', '001', 'ROI', 'converted', 'decon', 'ch01']
        
        if list0[-1]=="ch01":
            list1=list0[0:2]+["".join(list0[2]+list0[4])]+list0[5:]
            counter+=makeLink(newFolder, x, list1, ext)

# ex pour filename = 'scan_008_RT_9_10_001_ROI_converted_decon_ch02.tif'
# list0 = ['scan', '008', 'RT', '9', '10', '001', 'ROI', 'converted', 'decon', 'ch02']
# sinon si list0[-1]=="ch02"
# list1 = ['scan', '008', 'RT9', '001', 'ROI', 'converted', 'decon', 'ch02']
# list1 = [x if x!="ch02" else "ch01" for x in list1] remplace ch02 par ch01
# maintenant list1 = ['scan', '008', 'RT9', '001', 'ROI', 'converted', 'decon', 'ch01']

        elif list0[-1]=="ch02":
            list1=list0[0:2]+["".join(list0[2]+list0[3])]+list0[5:]
            list1=[x if x!="ch02" else "ch01" for x in list1]
            counter+=makeLink(newFolder, x, list1, ext)

# ex pour filename = 'scan_008_RT_9_10_001_ROI_converted_decon_ch00.tif'
# list0 = ['scan', '008', 'RT', '9', '10', '001', 'ROI', 'converted', 'decon', 'ch00']
# sinon si list0[-1]=="ch00"
# list1 = ['scan', '008', 'RT9', '001', 'ROI', 'converted', 'decon', 'ch00']

        elif list0[-1]=="ch00":
            list1=list0[0:2]+["".join(list0[2]+list0[3])]+list0[5:]
            counter+=makeLink(newFolder, x, list1, ext)

# ET list2 = ['scan', '008', 'RT10', '001', 'ROI', 'converted', 'decon', 'ch00']          
            list2=list0[0:2]+["".join(list0[2]+list0[4])]+list0[5:]
            counter+=makeLink(newFolder, x, list2, ext)
            
            
    print("{} RT links made for ROI {}".format(counter,ROI))