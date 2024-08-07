#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 13:13:10 2020

@author: christophe

"""
# Permet de copier des fichiers tiff d'une expérience (via un lien symbolique) 
# dans un répertoire choisi, pour ne pas toucher au dossier original.

import os, glob



rootFolder = "/mnt/grey/DATA/rawData_2020/Experiment_18_Olivier_Locus2L_low_res/Deconvolved/"
destFolder = "/home/christophe/Documents/Resultats/ROI_Olivier_Pourcentage_RT_By_Cell"

# Pour lister uniquement les dossiers présents dans le dossier rootFolder :
listFolder = [i for i in os.listdir(rootFolder) if os.path.isdir(rootFolder + os.sep + i)]
# 2eme possibilité (on obtient le chemin complet): 
#listFolder = glob.glob(rootFolder+ os.sep+ '*' + os.path.sep)


for folder in listFolder :
    filesRtDapi = glob.glob(rootFolder+ os.sep+ folder + os.sep + "*.tif")
    newFolder = destFolder + os.sep + folder
    os.mkdir(newFolder)
    for x in filesRtDapi :
        y = newFolder + os.sep + os.path.basename(x)
        os.system('ln -s ' + x + " " + y)
