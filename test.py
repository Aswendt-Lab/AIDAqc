# -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 15:37:07 2023

@author: arefks
"""
import os
import numpy as np
from nibabel.testing import data_path
import nibabel as nib
import matplotlib.pyplot as plt
from skimage import data
from skimage.filters import try_all_threshold
from skimage.filters import threshold_isodata
from scipy import ndimage
#%%
def snrCalclualtor2(input_file):

    #P = r"C:\Users\aswen\Documents\Data\Project_CRC\proc_data_processed\P1\StrokeBad\SP_T1_13_1_1_2_20180417_122304\T2w\SP_T1_13_1_1_2_20180417_122304.5.1.nii.gz"
    P = input_file
    #print(P)
    img = nib.load(P)
    Data = img.get_fdata()
    S = np.shape(np.squeeze(Data))
    if len(S) == 3:
        imgData = np.squeeze(Data)
    if len(S) == 4:
        imgData = np.squeeze(Data[:,:,:,int((S[2]/2))])
    
    
    #local thresholding
    imgData_new = np.zeros(S[0:3]);
    for ii in range(0,S[2]):
        temp_image = imgData[:,:,ii]
        global_thresh = threshold_isodata(temp_image)
        binary_global = temp_image > global_thresh
        imgData_new[:,:,ii] = binary_global
        
    
    COM=[int(i) for i in (ndimage.measurements.center_of_mass(imgData_new*imgData))]
    r = np.floor(0.10*(np.mean(S)))
    Mask = sphere(S, int(r) , COM)
    Singal = np.mean(imgData[Mask])
    
    
    x= int(S[0]*0.15)
    y = int(S[1]*0.15)
    z = int(S[2]*0.15)
    
    MaskN = np.zeros(S[0:3]);
    MaskN[:x,:y,:z] = 2
    MaskN[:x,-y:,:z] = 2
    MaskN[-x:,:y,:z] = 2
    MaskN[-x:,-y:,:z] = 2
    MaskN[:x,:y,-z:] = 2
    MaskN[:x,-y:,-z:] = 2
    MaskN[-x:,:y,-z:] = 2
    MaskN[-x:,-y:,-z:] = 2
    
    MaskND = nib.Nifti1Image(MaskN+Mask, a, img.header)
    nib.save(MaskND,os.path.join(r"C:\Users\aswen\Desktop\LabProject\TestSNR" ,os.path.basename(P) + '_Mask.nii.gz')) 
    nib.save(img, os.path.join(r"C:\Users\aswen\Desktop\LabProject\TestSNR" ,os.path.basename(P) + '_Data.nii.gz')) 
    
    
    n1 = np.squeeze(imgData[:x,:y,:z])
    n2 = np.squeeze(imgData[:x,-y:,:z])
    n3 = np.squeeze(imgData[-x:,:y,:z])
    n4 = np.squeeze(imgData[-x:,-y:,:z])
    n5 = np.squeeze(imgData[:x,:y,-z:])
    n6 = np.squeeze(imgData[:x,-y:,-z:])
    n7 = np.squeeze(imgData[-x:,:y,-z:])
    n8 = np.squeeze(imgData[-x:,-y:,-z:])
    
    
    Noise_std = np.std(np.array([n1,n2,n3,n4,n5,n6,n7,n8]))
    #show_slices([n8[:,:,3],np.squeeze(imgData[:,:,3])])
    #plt.show()
    SNR = Singal/Noise_std
    
    return SNR



def show_slices(slices):
   """ Function to display row of image slices """
   fig, axes = plt.subplots(1, len(slices))
   for i, slice in enumerate(slices):
       axes[i].imshow(slice.T, cmap="gray", origin="lower")
       

def sphere(shape, radius, position):
    """Generate an n-dimensional spherical mask."""
    # assume shape and position have the same length and contain ints
    # the units are pixels / voxels (px for short)
    # radius is a int or float in px
    assert len(position) == len(shape)
    n = len(shape)
    semisizes = (radius,) * len(shape)

    # genereate the grid for the support points
    # centered at the position indicated by position
    grid = [slice(-x0, dim - x0) for x0, dim in zip(position, shape)]
    position = np.ogrid[grid]
    # calculate the distance of all points from `position` center
    # scaled by the radius
    arr = np.zeros(shape, dtype=float)
    for x_i, semisize in zip(position, semisizes):
        # this can be generalized for exponent != 2
        # in which case `(x_i / semisize)`
        # would become `np.abs(x_i / semisize)`
        arr += (x_i / semisize) ** 2

    # the inner part of the sphere will have distance below or equal to 1
    return arr <= 1.0
#%%
SNR_list = []
for t in text_files:
    SNR = snrCalclualtor2(t)
    SNR_list.append(SNR)
    
df = pd.DataFrame()
df['FileAddress'] = text_files
df['SNR'] = SNR_list
df.to_csv(r"C:\Users\aswen\Desktop\LabProject\TestSNR\Result_SNR.csv")
