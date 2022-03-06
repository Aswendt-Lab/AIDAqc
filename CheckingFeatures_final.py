#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 25 13:00:02 2021

@author: kalantaria
"""

    #%% Loading nececcery libraries

import numpy as np
import QC
import pandas as pd
import pv_reader as pr
import openpyxl
import nibabel as nii
import os
import alive_progress as ap
     #%% Load Excel Data 

#%%
def main():
    import argparse
    
    parser = argparse.ArgumentParser(description='If you are about to use this tool, it means all raw data are already parsed.\
                                      Use this programm to calculate the features\
                                     in a seperate step as the second stage. The second stage checks\
                                     the features from every path found in the first stage according to its nature: for T2w and DTI:\
                                     SNR and resolution homogenity & for rsfmri tSNR, resoltion homogenity, Movement Severity based on mutual\
                                         information and Movement type. Movement type can be local or local & general. local movement happens ALWAYS.\
                                     if there is general movement detected, it indicates that the movement between the fmri images have a trend over time\
                                     , and are not just spontanious which means that the animal/patient/subject might has moved slowly in a specific di\
                                     rection. With kind regrads RFKS')
    parser.add_argument('address',help='Path of the excel file created in the first stage or any other sutibale path')
    args = parser.parse_args()
    Path = args.address    
    QC.CheckingFeatures(Path)
    
    
if __name__ == '__main__':
    main()
    
    #%% Plotting Result 
     
     # =============================================================================
     # 
     # f, axes = plt.subplots(2)
     # f.tight_layout(pad=3.0)
     # sns.distplot(snrBrum_vec, hist= True, ax= axes[0])
     # axes[0].set_title('Brum Method')
     # 
     # 
     # sns.distplot(snrCh_vec,hist = True, ax= axes[1])
     # axes[1].set_title('Chang Method')
     # 
     # for ax in axes.flat:
     #     ax.set(xlabel='SNR', ylabel='Quantity (%)')
     # 
     # =============================================================================
     