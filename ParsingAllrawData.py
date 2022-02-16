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
import CheckingFeatures_final as CF
#%% Command line interface
if __name__ == "__main__":
        
    parser = argparse.ArgumentParser(description='Parser of all MR files: Description:\
         This code will parse through every possible folder after a defined initial path,\
    looking for MR data files of any type. Then it will extract the wanted files \
    and eliminiate the duplicates.')
    parser.add_argument('initial_path', help='initial path to start the parsing (ending with "/" )')
    parser.add_argument('saving_path', help='Path where the results should be saved (ending with "/" )')
    parser.add_argument('-f','--forward',action='store_true',help='if this argument is used then the second stage of calculation will start right\
                        after the first stage automaticaly. The second stage checks\
                        the features from every path found in the first stage according to its nature: for T2w and DTI:\
                        SNR and resolution homogenity & for rsfmri tSNR, resoltion homogenity, Movement Severity based on mutual\
                            information and Movement type. Movement type can be local or local & general. local movement happens ALWAYS.\
                        if there is general movement detected, it indicates that the movement between the fmri images have a trend over time\
                        , and are not just spontanious which means that the animal/patient/subject might has moved slowly in a specific direction.\
                            \nWith kind regrads RFKS')                    
    parser.add_argument('-k', '--keyfolder', help='Definition keyfolder: A folder name which \
                        you are sure that every file has that folder in its path. You might ask: Why not state it \
                        in the initial path? The reason is that this is not always possible, see the following example:\
                        Initial path is "/Volumes/AG_Aswendt_Projects/ keyfolder is "MRI/" What the program does:\
                        /Volumes/AG_Aswendt_Projects/**/MRI/**/method" which means the program will start parsing every folder \
                        everything after AG_Aswendt_Progect AND it wont search for the method file unless it gets\
                        to a folder that is called MRI and after that it will parse for the file method. Note that \
                        if the keyfolder is not defined, the result will be the same, it just takes more TIME.\
                        ')                    
    args = parser.parse_args()
    initial_path = args.initial_path
    saving_path = args.saving_path
    keyfolder = args.keyfolder
    forward = args.forward
    #%% User Input: Main Path/folder where the program schould start to parse 
    
    #path = "/Volumes/AG_Aswendt_Projects/"
    Types = ['Dti*','EPI','RARE']
    
    
    #%% Parsing
    if not keyfolder:
        PathALL = initial_path + "**/method"
    else:
        PathALL = initial_path + "**/"+keyfolder+"**/method"
    
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
    
    Types_new = ['DTI','rsfMRI','T2w']
    for i,T in enumerate(Types_new):
        globals()['df'+ str(i)].to_excel(writer,sheet_name=T, index = False)
    
    
    dfError.to_excel(writer, sheet_name='ErrorData',index = False)
    
    writer.save()
    
    print('\n\nExcel file was created:' + str(saving_path2))
    print('\n\n%%%%%%%%%%%%%End of the first stage%%%%%%%%%%%%%%%'.upper())
    if forward == True:
        print('\nStarting Stage two ...'.upper())
        print('\nCalculating features...\n'.upper())
        print('This might take some time (hours/days) if the dataset is big enough!:) ...\n\n')
        CF.CheckingFeatures(saving_path2)
        