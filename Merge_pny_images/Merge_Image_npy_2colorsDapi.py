#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 10:59:21 2020

@author: christophe
"""


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 19 13:24:10 2020

@author: messina
"""
# project a number of images and makes an RGB image
import os
import glob
import matplotlib.pylab as plt
import numpy as np 
import cv2
import matplotlib.gridspec as gridspec
#from imageProcessing.imageProcessing import Image, imageAdjust
#import imageProcessing.makeProjections
from astropy.stats import SigmaClip
from photutils import Background2D, MedianBackground
from astropy.visualization import SqrtStretch, simple_norm
import operator
from PIL import Image
from PIL import ImageDraw
import numpy
from PIL import Image,ImageEnhance
#%% loads RTs images from alignImages folder
rootFolderRT = "/mnt/grey/DATA/rawData_2021/Experiment_13_Tof/Images_deconvoled/ROI009/alignImages"
fileNameRTa = 'scan_002_RT1_009_ROI_converted_decon_ch01_2d_registered.npy'
fileNameRTb = 'scan_002_RT4_009_ROI_converted_decon_ch01_2d_registered.npy'

filename_RT_a = rootFolderRT + os.sep + fileNameRTa # Le RT_a sera en rouge
filename_RT_b = rootFolderRT + os.sep + fileNameRTb # Le RT_b sera en vert

# loads Dapi image from segmentedObjects folder
rootFolderDapiMask="/mnt/grey/DATA/rawData_2021/Experiment_13_Tof/Images_deconvoled/ROI009/segmentedObjects"
fileNameDapiMask = 'scan_001_DAPI_009_ROI_converted_decon_ch00_Masks.npy'
filename_mask = rootFolderDapiMask + os.sep + fileNameDapiMask


Barcode_a = np.load(filename_RT_a)
Barcode_b = np.load(filename_RT_b)
Mask_Dapi = np.load(filename_mask)
#Fluoreseine = Image.open(filename_fluoresceine)
#Mask_normalized = Mask/Fluoreseine

rgbArray = np.zeros((2048,2048,3), 'uint8') # mettre 2 ou 3 si 2 ou 3 images a superposer
rgbArray[..., 0] = (Barcode_a)
rgbArray[..., 1] = (Barcode_b)
rgbArray[..., 2] = (Mask_Dapi)
img = Image.fromarray(rgbArray)

# substract background image BarcodeA
sigma_clip = SigmaClip(sigma=3.0)
bkg_estimator = MedianBackground()
bkg = Background2D(Barcode_a,(64,64),filter_size=(3,3),sigma_clip=sigma_clip,bkg_estimator=bkg_estimator)
Barcode_a_substracted = Barcode_a - bkg.background
Barcode_a_substracted = Barcode_a_substracted - np.min(Barcode_a_substracted)

# # substract background image BarcodeB
sigma_clip = SigmaClip(sigma=3.0)
bkg_estimator = MedianBackground()
bkg = Background2D(Barcode_b,(64,64),filter_size=(3,3),sigma_clip=sigma_clip,bkg_estimator=bkg_estimator)
Barcode_b_substracted = Barcode_b - bkg.background
Barcode_b_substracted = Barcode_b_substracted - np.min(Barcode_b_substracted)

# substract background image Mask
sigma_clip = SigmaClip(sigma=3.0)
bkg_estimator = MedianBackground()
bkg = Background2D(Mask_Dapi,(64,64),filter_size=(3,3),sigma_clip=sigma_clip,bkg_estimator=bkg_estimator)
Mask_Dapi_substracted = Mask_Dapi - bkg.background
Mask_Dapi_substracted = Mask_Dapi_substracted - np.min(Mask_Dapi_substracted)

# Test normalization Mask 

Barcode_a_substracted = ((Barcode_a_substracted-Barcode_a_substracted.min())/(Barcode_a_substracted.max()-Barcode_a_substracted.min()))
Barcode_b_substracted = ((Barcode_b_substracted-Barcode_b_substracted.min())/(Barcode_b_substracted.max()-Barcode_b_substracted.min()))
Mask_Dapi_substracted = ((Mask_Dapi_substracted-Mask_Dapi_substracted.min())/(Mask_Dapi_substracted.max()-Mask_Dapi_substracted.min()))

# ATTENTION la matrice normalisée a des valeurs comprises entre 0 et 1.
# La matrice RGB a un max d'intensité de 255. Il faudra donc multiplier au maximum la valeur 
# du pixel Barcode (qui est comprise entre 0 et 1) par 255, sinon les valeurs max de 1 seront abérantes
rgbArray = np.zeros((2048,2048,3), 'uint8')
rgbArray[..., 0] = (Barcode_a_substracted)*500 #255 max (rouge)
rgbArray[..., 1] = (Barcode_b_substracted)*800 #255 max (vert)
rgbArray[..., 2] = (Mask_Dapi_substracted)*200 #255 max

#min = 0, max = 2048 
img = Image.fromarray(rgbArray)
plt.imshow(img)
plt.xlim(0,2048)
plt.ylim(0,2048)
plt.show()

# img = Image.fromarray(rgbArray)
# plt.imshow(img)
# plt.xlim(860,905)
# plt.ylim(1135,1178)
# plt.show()