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
#%% Command line interface
if __name__ == "__main__":
        
    parser = argparse.ArgumentParser(description='Parser of all MR files: Description:\
         This code will parse through every possible folder after a defined initial path,\
    looking for MR data files of any type. Then it will extract the wanted files \
    and eliminate the duplicates.')
    parser.add_argument('initial_path', help='initial path to start the parsing (ending with "/" )')
    parser.add_argument('saving_path', help='Set the path where the results should be saved (ending with "/" )')
    parser.add_argument('-f','--forward',action='store_true',help='Set this parameter if you JUST want to parse the databank without the\
                        the quality measurement (Just creates an excel table with all existing MR files addresses relating to DTI, T2w and fMRI sequences)')                                       
    parser.add_argument('-e','--exclude',type=str, choices=['T2w', 'fMRI', 'DTI'],help='Set this parameter if you want to \
                        exclude a specific type of sequence between the three types of T2w, fMRI and DTI. Example use \
                            in terminal: python ParsingAllrawData.py <inital_path>/ <saving_path>/ -e fMRI')
    args = parser.parse_args()
    initial_path = args.initial_path
    saving_path = args.saving_path
    forward = args.forward
    exclude = args.exclude
    #%% User Input: Main Path/folder where the program schould start to parse 
    #path = "/Volumes/AG_Aswendt_Projects/"
    
    Types = ['Dti*','EPI','RARE']
    Types_new = ['DTI','rsfMRI','T2w']
    
    if exclude == 'DTI':
        Types.remove('Dti*')
        Types_new.remove('DTI')
        
    if exclude == 'fMRI':
        Types.remove('EPI')
        Types_new.remove('rsfMRI')
        
    if exclude == 'T2w':
        Types.remove('RARE')
        Types_new.remove('T2w')
    QC.tic()
    print("Hello!")
    print('------------------------------------------------------------')
    print('Thank you for using our Code. For questions please contact us over:')
    print('aref.kalantari-sarcheshmeh@uk-koeln.de or markus.aswendt@uk-koeln.de')
    print('Lab: AG Neuroimaging and neuroengineering of experimental stroke University Hospital Cologne')
    print('Web: https://neurologie.uk-koeln.de/forschung/ag-neuroimaging-neuroengineering/')
    print('------------------------------------------------------------')
        
    
    #%% Parsing
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
    
    #%% Saving parsed files to excel sheets
    
    
    for i,T in enumerate(Types):
        globals()['df'+ str(i)] = pd.DataFrame(ABook[T])
    
    
    dfError = pd.DataFrame()
    dfError['ErrorData'] = ErrorList
    
    #saving_path_temp = '/Volumes/AG_Aswendt_Share/Scratch/Aref/Results/AGAswendtMRData.xlsx'
    saving_path2 = saving_path + 'QuiC_Data_Result.xlsx'
    writer = pd.ExcelWriter(saving_path2, engine='xlsxwriter')
    
    #ABook.keys()
    
    for i,T in enumerate(Types_new):
        globals()['df'+ str(i)].to_excel(writer,sheet_name=T, index = False)
    
    
    dfError.to_excel(writer, sheet_name='ErrorData',index = False)
    
    writer.save()
    
    print('\n\nExcel file was created:' + str(saving_path2))
    print('\n\n%%%%%%%%%%%%%End of the first stage%%%%%%%%%%%%%%%'.upper())
    if forward == True:
        print("***")
    else:
        print('\nStarting Stage two ...'.upper())
        print('\nCalculating features...\n'.upper())
        print('This might take some time (hours/days) if the dataset is big enough!:) ...\n\n')
        QC.CheckingFeatures(saving_path2)
        QC.toc()
        print('------------------------------------------------------------')
        print('Thank you for using our Code. For questions please contact us over:')
        print('aref.kalantari-sarcheshmeh@uk-koeln.de or markus.aswendt@uk-koeln.de')
        print('Lab: AG Neuroimaging and neuroengineering of experimental stroke University Hospital Cologne')
        print('Web:https://neurologie.uk-koeln.de/forschung/ag-neuroimaging-neuroengineering/')
        print('------------------------------------------------------------')
        