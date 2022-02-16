#%% 
"""
Name: Aref Kalantari
Email: arefks@gmail.com
Date: 24.08.2021
-----------------------------
Code Describtion: Quality Control. This code calculates various factors and parameters of MRI images, combining them to
a final grade form 0 to 20, with 0 standting for Worst Quality and 20 for the best Quality. For autmation porpuses the 
idea is to use datasets with a grade higiher then 12 and to discard the rest. For semiautomatic purposes the resulting 
parameters are all observable to the user and the user can decide to continue with the dataset or to dicard it.
-----------------------------

"""

#%% Loading nececcery libraries

import numpy as np
import nibabel as nii
import changSNR as ch
import brummerSNR as bm
import glob
import MI as mi
import BrukerMRI as br
from matplotlib.pyplot import imshow

#%% SNR functions

def ResCalculator(input_file):
    
    HDR = input_file.header
    Spati_Res = HDR['pixdim'][1:4]
    
    return Spati_Res

def snrCalclualtor(input_file):

    imgData = input_file
    ns = imgData.shape[2]  # Number of slices
    nd = imgData.ndim
    ns_lower = int(np.floor(ns/2) - 2)
    ns_upper = int(np.floor(ns/2) + 2)
    
    noiseChSNR = np.zeros(ns)
    IM = np.asanyarray(imgData.dataobj)
    imgData = np.ndarray.astype(IM, 'float64')
    #print('/NewData/',end=" ")
    for slc in range(ns_lower,ns_upper):
        #   Print % of progress
        #print('S' + str(slc + 1), end=",")

        # Decision if the input data is DTI type or T2w
        if nd == 3:
            slice = imgData[:, :, slc]
        
        if nd == 4:
            slice = imgData[:, :, slc,1]
            
        curSnrCHMap, estStdChang, estStdChangNorm = ch.calcSNR(slice, 0, 1)
        noiseChSNR[slc] = estStdChang

    snrCh = 20 * np.log10(np.mean(imgData) / np.mean(noiseChSNR))

    return snrCh



def TsnrCalclualtor(input_file):
    
    imgData = input_file
    IM = np.asanyarray(imgData.dataobj)
    imgData = np.ndarray.astype(IM, 'float64')
    temp_mean = imgData.mean(axis=(0,1,3))
    temp_max = temp_mean.argmax()

    ns = imgData.shape[2]  # Number of slices
    nt = imgData.shape[-1]
    nd = imgData.ndim
    
    if temp_max == 0:
        temp_max = temp_max+1
        print('temp_max is the min slice')
        
    if temp_max == ns:
        temp_max = temp_max-1
        print('temp_max is the max slice')
        
        
    ns_lower = int(temp_max - 1)
    ns_upper = int(temp_max + 1)

    noiseChSNR = np.zeros(ns)
    tSNR_Ch_vec = []
    for slc in range(ns_lower,ns_upper):
        snrCh_t = []
        
        for t in range(1,nt):
            
            slice = imgData[:, :, slc,t]
            
            curSnrCHMap, estStdChang, estStdChangNorm = ch.calcSNR(slice, 0, 1)
                
            noiseChSNR = estStdChang
            snrCh_temp = 20 * np.log10(np.mean(imgData[:,:,slc,t]) / noiseChSNR)
            snrCh_t.append(snrCh_temp)
    
        Std_Ch = np.std(snrCh_t)
        tSNR_Ch = np.mean(snrCh_t)/Std_Ch
        tSNR_Ch_vec.append(tSNR_Ch)
        
        
    tSNR_Ch_final = np.mean(tSNR_Ch_vec)
   
    return tSNR_Ch_final



def Ismovement(input_file):
    TypeMov=[]
    imgData = input_file
    IM = np.asanyarray(imgData.dataobj)
    imgData = np.ndarray.astype(IM, 'float64')
    S = np.shape(imgData)
    temp_mean = imgData.mean(axis=(0,1,3))
    temp_max = temp_mean.argmax()
    temp_Data = imgData[:,:,temp_max,:]
    Im_fix = temp_Data[:,:,0]
    Im_rot = temp_Data
    
    MI_all = []
    for z in range(S[-1]):
        
        MI = mi.mutualInfo(Im_fix,Im_rot[:,:,z])
        MI_all.append(MI)
    
    Final = np.asarray(MI_all)
    
    m,b = np.polyfit(np.arange(0,len(Final)),Final,1)
    if m > 0.2:
        TypeMov = 'General & Local'
    else:
        TypeMov = 'Local'
        
    GMV = getrange(Final)
    LMV = 3*np.std(Final)
    
    return Final,GMV,LMV,TypeMov

def getrange(numbers):
    return max(numbers) - min(numbers)






