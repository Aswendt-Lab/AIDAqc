#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 25 13:00:02 2021

@author: kalantaria
"""

    #%% Loading nececcery libraries

import numpy as np
import QC
import time
import pandas as pd
import pv_reader as pr
import openpyxl
from matplotlib.pyplot import imshow
import nibabel as nii
import os
import alive_progress as ap


     #%% Load Excel Data 
def CheckingFeatures(Path):   
    
  
    Names = []
    ErorrList = []
    #Path = "/Users/kalantaria/Desktop/Res/QuiC_Data_Result2.xlsx" #!!!
    xls = pd.ExcelFile(Path,engine= 'openpyxl')
    Names = xls.sheet_names
    Names_temp = Names
    saving_path = os.path.dirname(Path) 
    
    if 'ErrorData' in Names:
        Names.remove('ErrorData')
    
    Abook = []
    for n in Names:
        Abook.append(pd.read_excel(Path,engine= 'openpyxl',sheet_name = n))
    
    Names = Names_temp
    #%% Calculate all the SNR for all the data that was found 
    # in the last step and saving it into a vector
    # Load Bruker data from each address 
    saving_path2 = saving_path + '/QuiC_Data_Result_Processed_featurs_testing.xlsx'
    writer = pd.ExcelWriter(saving_path2, engine='xlsxwriter')
    kk =0
    for ii,N in enumerate(Names):
        if N != 'ErrorData':
            if kk > 0:
                print(str(kk) + 'faulty files were found:All faulty files are available in the\
                      Errorlist tab in the Excel outputs\n')
            print(N+' processing... \n')
            text_files = Abook[ii][0]
            snrCh_vec =[]
            SpatRes_vec = []
            MI_vec_all = []
            LMV_all = []
            TypeMov_all = []
            text_files_new = []
            kk = 0
            i=1
            
            with ap.alive_bar(len(text_files),spinner='wait',refresh_secs = 0) as bar:
                for tf in text_files:
                    
                    path_split = tf.split('/')
                    
                    procno = str(1)
                    expno = path_split[-1]
                    study = path_split[-2]
                    raw_folder = '/'.join(path_split[:-2])
                    proc_folder = raw_folder+ '/proc_data' #Here still adjustment is needed
                    pv = pr.ParaVision(proc_folder, raw_folder, study, expno, procno)
                    
                    CP_v = tf + '/pdata/1/visu_pars' # Check Parameter: Visu_pars
                    CP_a = tf + '/acqp' # Check Parameter: acqp
                    
                    if os.path.isfile(CP_v) and os.path.isfile(CP_a):
                        pv.read_2dseq( map_raw=False, map_pv6=False, roll_fg=False, squeeze=False, compact=False, swap_vd=False, scale=1.0)
                        input_file = nii.squeeze_image(pv.nifti_image)
                    else:
                        ErorrList.append(tf)
                        continue
                   
                    # Resoultution
                    SpatRes = QC.ResCalculator(input_file)
                    
                    
                    if N == 'T2w' or N == 'DTI':
                        # Signal 2 noise ratio
                        try:
                            snrCh = QC.snrCalclualtor(input_file)
                        except Exception:
                            ErorrList.append(tf)
                            kk = kk+1
                            continue
                            
                        LMV_all = np.nan
                        TypeMov_all = np.nan
                        
                    if N == 'rsfMRI':
                        #temporal signal 2 noise ratio
                        try:
                            snrCh = QC.TsnrCalclualtor(input_file)
                        except Exception:
                            ErorrList.append(tf)
                            kk = kk+1
                            continue
                        
                        # movement severity with the help of mutual information
                        Final,GMV,LMV,TypeMov = QC.Ismovement(input_file)
                        TypeMov_all.append(TypeMov)
                        LMV_all.append(LMV)
                        MI_vec_all.append(Final)
                        
                    SpatRes_vec.append(SpatRes)      
                    snrCh_vec.append(snrCh)
                    i=i+1
                    text_files_new.append(tf)
                    bar()
                    
            # Saving parsed files to excel sheets
            AR = [text_files_new,np.array(SpatRes_vec),np.array(snrCh_vec),np.array(LMV_all),TypeMov_all]
            
            # using the savetxt 
            # from the numpy module
            
            df = pd.DataFrame()
            df['FileAddress'] = AR[0]
            df['SpatRx'] = AR[1][:,0]
            df['SpatRy'] = AR[1][:,1]
            df['Slicethick'] = AR[1][:,2]
            
            if N == 'T2w' or N == 'DTI':
                 df['SNR Chang'] = AR[2]
                 
            else:
                 df['tSNR Chang'] = AR[2]
                 df['Local Movement Variability']=AR[3]
                 df['Movement Type']=AR[4]
                 
            
            
             
            df.to_excel(writer,sheet_name=N, index = False)
        
        else:
            df = pd.DataFrame()
            df['ErorrList'] = ErorrList
            df.to_excel(writer,sheet_name=N, index = False)
        
    writer.save()
    print('\n\nExcel file was created:' + str(saving_path2))
    
    print('\n\n%%%%%%%%%%%%%End of the Second stage%%%%%%%%%%%%%%%\n\n'.upper())
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
    CheckingFeatures(Path)
    
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
     