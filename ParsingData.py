#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 22 17:28:14 2021
@author: kalantaria
Description: This code will parse through every possible folder after a defined initial path,
looking for MR data files of any type. Then it will extract the wanted files 
and eliminiate the duplicates.
"""

import os
import glob
import pv_parser as par
import re
import pandas as pd
import argparse
import alive_progress as ap
import numpy as np
import QC
import time
from openpyxl import Workbook
#%% Command line interface
if __name__ == "__main__":
        
    parser = argparse.ArgumentParser(description='Parser of all MR files: Description:\
          This code will parse through every possible folder after a defined initial path,\
     looking for MR data files of any type. Then it will extract the wanted files \
     and eliminate any duplicates.')
    parser.add_argument('-i','--initial_path',required=True \
                        help='initial path to start the parsing (ending with "/" )')
    parser.add_argument('-o','--output_path',required=True,\
                        help='Set the path where the results should be saved (ending with "/" )')
    parser.add_argument('-f','--format_type',\
                        help="you need to tell what kind of format your images are :\
                            nifti or raw",type=str,required=True,choices=["nifti","raw"])  
    parser.add_argument('-t','--sequence_types',\
                        help="you need to tell what kind of Sequences should be used in \
                             for processing the dataset:\
                            T2w, DTI, fmri",type=str,required=False,choices=["T2w","DTI","fMRI"],default=["T2w","DTI","fMRI"])  
    parser.add_argument('-s','--suffix',\
                        help="If necessery you can specify what kind of sufix the data to look for should have :\
                            for example: -s test , this means it will only look for data that have this\
                                suffix befor the .nii.gz, meaing test.nii.gz",type=str,required=False, Default=None)  
                                                  
    
    
    
    
    
    args = parser.parse_args()
    initial_path = args.initial_path
    saving_path = args.output_path
    format_type= args.format_type
    sequence_types = args.sequence_types
    suffix = args.suffix
    #%% User Input: Main Path/folder where the program schould start to parse 
    #path = "/Volumes/AG_Aswendt_Projects/"
    
    Types = ['Dti*','EPI','RARE']
    Types_new = ['DTI','rsfMRI','T2w']
    
    # if exclude == 'DTI':
    #     Types.remove('Dti*')
    #     Types_new.remove('DTI')
        
    # if exclude == 'fMRI':
    #     Types.remove('EPI')
    #     Types_new.remove('rsfMRI')
        
    # if exclude == 'T2w':
    #     Types.remove('RARE')
    #     Types_new.remove('T2w')
    QC.tic()
    print("Hello!")
    print('------------------------------------------------------------')
    print('Thank you for using our Code. For questions please contact us over:')
    print('aref.kalantari-sarcheshmeh@uk-koeln.de or markus.aswendt@uk-koeln.de')
    print('Lab: AG Neuroimaging and neuroengineering of experimental stroke University Hospital Cologne')
    print('Web: https://neurologie.uk-koeln.de/forschung/ag-neuroimaging-neuroengineering/')
    print('------------------------------------------------------------')

    #%% Parsing
    if format_type== "raw":
        PathALL = initial_path + "**/method"
        with ap.alive_bar(title='Parsing through folders ...',length=10,stats = False,monitor=False) as bar:
            text_files = glob.glob(PathALL, recursive = True)
            kall = len(text_files)
        
        print(( 'Total number of '+ str(kall) + ' files were found:'+'Parsing finished! '.upper()).upper())
    
        #%% Extrtacting usable data
        ABook = {}
        ErrorList =[]
        CheckDates = []
        C = 0
        for i in range(len(Types)):         #Creation of Adress Book
            ABook[Types[i]] = []
    
    
        with ap.alive_bar(kall, title='Extracting T2w, DTI and fmri files:'.upper(),length=10,stats = False,spinner= 'wait') as bar:   
            for p in text_files:   #filling the Address Book with wanted files
            
                try:
                
                    NameTemp = par.read_param_file(p)
                    MN = NameTemp[1]["Method"]
                    DateTemp = NameTemp[0]['Date']
                
                #  NameTemp = par.read_param_file(p+'/acqp')
                #  MN2 = NameTemp[1]["ACQ_scan_name"]
                
                    Ans = []
                except SystemExit:
                    ErrorList.append(p)
                
            
                if DateTemp not in CheckDates:
                
                    for i,t in enumerate(Types):
                    
                        typ = re.search(t,MN)
                        if typ != None:
                            Ans = i
                
                    if Ans != []:
                        ABook[Types[Ans]].append(os.path.dirname(p))
                        C = C+1
                    
                CheckDates.append(DateTemp)
                bar()
        M = dict.fromkeys(CheckDates)
        
        print(' '+str(C)+' files were extracted! %%%'.upper())
        print((' ' + str(len(CheckDates)-len(M))+ ' Duplicates were Eliminated! %%%').upper())
        print("a",ABook)
        #%% Saving parsed files 
       
        #saving in csv file






        dti_addreses= pd.DataFrame(ABook["Dti*"])
        dti_path= os.path.join(saving_path,"raw_data_addreses_DTI.csv")
        dti_addreses.to_csv(dti_path, sep=',',index=False)

        t2w_addreses= pd.DataFrame(ABook["RARE"])
        t2w_path= os.path.join(saving_path,"raw_data_addreses_T2w.csv")
        t2w_addreses.to_csv(t2w_path, sep=',',index=False)

        rsfmri_addreses= pd.DataFrame(ABook["EPI"])
        rsfmri_path= os.path.join(saving_path,"raw_data_addreses_fMRI.csv")
        rsfmri_addreses.to_csv(rsfmri_path, sep=',',index=False)


        print('\n\ncsv files were created:' + str(saving_path))
        print('\n\n%%%%%%%%%%%%%End of the first stage%%%%%%%%%%%%%%%'.upper())









        # to make exel file as an output you can uncomment below lines
        # for i,T in enumerate(Types):
        #     globals()['df'+ str(i)] = pd.DataFrame(ABook[T])
    
    
        # dfError = pd.DataFrame()
        # dfError['ErrorData'] = ErrorList
    
       
        # saving_path2 = saving_path + 'QuiC_Data_Result_raw.xlsx'
        # writer = pd.ExcelWriter(saving_path2, engine='xlsxwriter')
    
        
    
        # for i,T in enumerate(Types_new):
        #     globals()['df'+ str(i)].to_excel(writer,sheet_name=T, index = False)
    
    
        # dfError.to_excel(writer, sheet_name='ErrorData',index = False)
    
        # writer.save()
    
        # print('\n\nExcel file was created:' + str(saving_path2))
        # print('\n\n%%%%%%%%%%%%%End of the first stage%%%%%%%%%%%%%%%'.upper())        
    
    #%% Parsing nifti format

    if format_type=="nifti":

        PathALL = initial_path + os.path.join("**","*1.nii.gz")
        with ap.alive_bar(title='Parsing through folders ...',length=10,stats = False,monitor=False) as bar:
            text_files = glob.glob(PathALL, recursive = True)
            kall = len(text_files)
        
        print(( 'Total number of '+ str(kall) + ' files were found:'+'Parsing finished! '.upper()).upper())
        ABook={}
        for i in range(len(Types)):         #Creation of Adress Book
            ABook[Types_new[i]] = []
    

        for i,T in enumerate(Types_new):
            globals()['df'+ str(i)] = pd.DataFrame(ABook[T])        
        for i in text_files :

            if "DTI"in i :
                ABook["DTI"].append(i)
            elif "rsfMRI"  in i :
                ABook["rsfMRI"].append(i)
            elif "fMRI" in i :
                ABook["rsfMRI"].append(i)
            elif "T2w"in i :
                ABook["T2w"].append(i)

        #saving in csv file
        dti_addreses= pd.DataFrame(ABook["DTI"])
        dti_path= os.path.join(saving_path,"nifti_data_addreses_DTI.csv")
        dti_addreses.to_csv(dti_path, sep=',',index=False)

        t2w_addreses= pd.DataFrame(ABook["T2w"])
        t2w_path= os.path.join(saving_path,"nifti_data_addreses_T2w.csv")
        t2w_addreses.to_csv(t2w_path, sep=',',index=False)

        rsfmri_addreses= pd.DataFrame(ABook["rsfMRI"])
        rsfmri_path= os.path.join(saving_path,"nifti_data_addreses_fMRI.csv")
        rsfmri_addreses.to_csv(rsfmri_path, sep=',',index=False)


        print('\n\ncsv files were created:' + str(saving_path))
        print('\n\n%%%%%%%%%%%%%End of the first stage%%%%%%%%%%%%%%%'.upper())
       
       
        
        # to make exel file as an output you can uncomment below lines


        # for i,T in enumerate(Types_new):
        #     globals()['df'+ str(i)] = pd.DataFrame(ABook[T])
    
    

    
        
        saving_path2 = saving_path + 'QuiC_Data_Result_nifti.xlsx'
        # writer = pd.ExcelWriter(saving_path2, engine='xlsxwriter')
    
       
    
        # for i,T in enumerate(Types_new):
        #     globals()['df'+ str(i)].to_excel(writer,sheet_name=T, index = False)
    
    
        
    
        # writer.save()
        
        #print('\n\nExcel file was created:' + str(saving_path2))


        
    if forward == True:
        print("***")
    else:
        print('\nStarting Stage two ...'.upper())
        print('\nCalculating features...\n'.upper())
        print('This might take some time (hours/days) if the dataset is big enough!:) ...\n\n')
        if format_type=="raw":

            QC.CheckingrawFeatures(saving_path,dti_addreses,rsfmri_addreses,t2w_addreses)
            QC.toc()
        elif format_type=="nifti":
            QC.CheckingNiftiFeatures(saving_path,dti_path,rsfmri_path,t2w_path)
            QC.toc()

            


        print('------------------------------------------------------------')
        print('Thank you for using our Code. For questions please contact us over:')
        print('aref.kalantari-sarcheshmeh@uk-koeln.de or markus.aswendt@uk-koeln.de')
        print('Lab: AG Neuroimaging and neuroengineering of experimental stroke University Hospital Cologne')
        print('Web:https://neurologie.uk-koeln.de/forschung/ag-neuroimaging-neuroengineering/')
        print('------------------------------------------------------------')








    
    

    
    
    
    
    
    
    
    
    