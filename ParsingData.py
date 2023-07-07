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
import pandas as pd
import argparse
import alive_progress as ap
import QC
import shutil
#from openpyxl import Workbook
import FeatureCheck as fc
#%% Command line interface
if __name__ == "__main__":
        
    parser = argparse.ArgumentParser(description='Parser of all MR files: Description:\
          This code will parse through every possible folder behind a defined initial path,\
     looking for MR data files of any type. Then it will extract the wanted files \
     and eliminate any duplicates(ex:python ParsingData.py -i C:\BME\aida\raw_data -o C:\BME\aida\raw_data -f raw.')
    parser.add_argument('-i','--initial_path',required=True, \
                        help='initial path to start the parsing')
    parser.add_argument('-o','--output_path',required=True,\
                        help='Set the path where the results should be saved')
    parser.add_argument('-f','--format_type',\
                        help="the format your dataset has:\
                            nifti or raw Bruker",type=str,required=True,choices=["nifti","raw"])  
   # parser.add_argument('-t','--sequence_types',\
   #                     help="you need to tell what kind of Sequences should be used in \
   #                          for processing the dataset:\
   #                         T2w, DTI, fmri",type=str,required=False,choices=["T2w","DTI","fMRI"],default=["T2w","DTI","fMRI"])  
    parser.add_argument('-s','--suffix',\
                        help="If necessary you can specify what kind of suffix the data to look for should have :\
                            for example: -s test , this means it will only look for data that have this\
                                suffix befor the .nii.gz, meaning test.nii.gz",type=str, required=False, default="")  
                                                  
    
    
    args = parser.parse_args()
    initial_path = args.initial_path
    saving_path = args.output_path
    format_type= args.format_type
    #sequence_types = args.sequence_types
    suffix = args.suffix
    #%% Information for the user 
    
    QC.tic()
    print("Hello! Are you ready to get rid of bad quality data?")
    print('------------------------------------------------------------')
    print('Thank you for using our code. Contact:')
    print('aref.kalantari-sarcheshmeh@uk-koeln.de / markus.aswendt@uk-koeln.de')
    print('Lab: AG Neuroimaging and Neuroengineering of Experimental Stroke, University Hospital Cologne')
    print('Web: https://neurologie.uk-koeln.de/forschung/ag-neuroimaging-neuroengineering/')
    print('------------------------------------------------------------')
    #%% Path Construction
    
    if not os.path.exists(saving_path):
        os.mkdir(saving_path)
    
    #%% Parsing
    
    Types = ['Dti','EPI','RARE']
    Types_new = ['DTI','rsfMRI','T2w']
  
    if format_type == "raw":
        
        
        
        PathALL = os.path.join(initial_path,"**","acqp")
        with ap.alive_bar(title='Parsing through folders ...',length=10,stats = False,monitor=False) as bar:
            text_files = glob.glob(PathALL, recursive = True)
            kall = len(text_files)
        
        print(( 'Total number of '+ str(kall) + ' files were found:'+' parsing finished! '.upper()).upper())
    
        #%% Extrtacting usable data
        ABook = {}
        ErrorList =[]
        CheckDates = []
        C = 0
        DTI_string = ["DTI","STRUCT","DWI"]
        FMRI_string = ["RESTING","FUN","RS","FMRI","BOLD"]
        T2_string = ["T2","T1","ANAT","RARE","TURBO"]
        NotAllowed = ["LOC","PIL","FISP","MAP","WOB","NOIS"]
        #EPI_flag = ["EPI"]
        
        
        for i in range(len(Types)):         #Creation of Adress Book
            ABook[Types[i]] = []
    
    
        with ap.alive_bar(kall, title='Extracting T1 or T2 weighted, structural and functional sequences:'.upper(),length=10,stats = False,spinner= 'wait') as bar:   
            
            for p in text_files:   #filling the Address Book with wanted files
            
                try:
                
                    NameTemp = par.read_param_file(p)
                    MN = NameTemp[1]["ACQ_method"].upper()  #Here we check what the name of the sequence is
                    MN2 = NameTemp[1]["ACQ_protocol_name"].upper()
                    KEY = MN + MN2
                    DateTemp = NameTemp[0]['Date'] #Here we check the date of the measurement
                    Ans = []
                except KeyError:
                    print(p)
                    ErrorList.append(p)
                
                Flag_anat = any([(aa in KEY) for aa in T2_string])
                Flag_struct = any([(aa in KEY) for aa in DTI_string])
                Flag_func = any([(aa in KEY) for aa in FMRI_string])
                Flag_notAllowed = any([(aa in KEY) for aa in NotAllowed])
                Flag_epi = "EPI" in KEY
                
            
                if DateTemp not in CheckDates:
                    
                    if Flag_struct and not Flag_notAllowed:
                        ABook["Dti"].append(os.path.dirname(p))
                        C = C+1
                    elif Flag_func and not Flag_notAllowed:
                        ABook["EPI"].append(os.path.dirname(p)) #I know it is totally confusing with EPI as the col name for the ABook but sadly EPI can also be a DTI scan
                        C = C+1
                    elif Flag_anat and not Flag_notAllowed:
                        ABook["RARE"].append(os.path.dirname(p))
                        C = C+1
                    elif Flag_epi and not Flag_notAllowed:
                        TP = NameTemp[1]["ACQ_time_points"]
                        if max(TP) == len(TP)-1 and any(TP):
                            ABook["EPI"].append(os.path.dirname(p))
                            C = C+1
                        elif any(TP):
                            ABook["Dti"].append(os.path.dirname(p)) #I know it is totally confusing with EPI as the col name for the ABook but sadly EPI can also be a DTI scan
                            C = C+1
                        
                        
 #                   for i,t in enumerate(Types):
 #                   
                       # if t in MN or t in p:
 #                        if t.upper() in MN.upper():
 #                           ABook[Types[i]].append(os.path.dirname(p))
 #                          C = C+1
                    
                CheckDates.append(DateTemp)
                bar()
                
        M = dict.fromkeys(CheckDates)
        
        print(' '+str(C)+' files were extracted! %%%'.upper())
        print((' ' + str(len(CheckDates)-len(M))+ ' duplicates were eliminated! %%%').upper())
        #%% Saving parsed files 
       
        #saving in csv file
        for n,type in enumerate(Types):
             if len(ABook[type]) !=0:
                 addreses= pd.DataFrame(ABook[type])
                 csv_path= "raw_data_addreses_"+Types_new[n]+".csv"
                 csv_path= os.path.join(saving_path,csv_path)
                 addreses.to_csv(csv_path, sep=',',index=False)
        if ErrorList:
            dfError = pd.DataFrame()
            dfError['ErrorData'] = ErrorList
            eror= os.path.join(saving_path,"CanNotOpenTheseFiles.csv")
            dfError.to_csv(eror,index=False)

        print('\n\ncsv files were created:' + str(saving_path))
        print('\n\n%%%%%%%%%%%%% End of stage 1 %%%%%%%%%%%%%%%'.upper())

        # to make exel file as an output you can uncomment below lines
        # for i,T in enumerate(Types):
        #     globals()['df'+ str(i)] = pd.DataFrame(ABook[T])
    
    

    
       
        # saving_path2 = saving_path + 'QuiC_Data_Result_raw.xlsx'
        # writer = pd.ExcelWriter(saving_path2, engine='xlsxwriter')
    
        
    
        # for i,T in enumerate(Types_new):
        #     globals()['df'+ str(i)].to_excel(writer,sheet_name=T, index = False)
    
    
        # dfError.to_excel(writer, sheet_name='ErrorData',index = False)
    
        # writer.save()
    
        # print('\n\nExcel file was created:' + str(saving_path2))
        # print('\n\n%%%%%%%%%%%%%End of the first stage%%%%%%%%%%%%%%%'.upper())        
    
    #%% Parsing nifti format

    elif format_type=="nifti":

        DTI_string = ["DTI","STRUCT","DWI"]
        FMRI_string = ["RESTING","FUN","RS","FMRI","BOLD"]
        T2_string = ["T2W","T1W","ANAT","RARE","TURBO"]
        NotAllowed = ["LOC","PIL","FISP","MAP","WOB"]

        PathALL = os.path.join(initial_path,"**","*" + suffix + ".nii*")
        with ap.alive_bar(title='Parsing through folders ...',length=10,stats = False,monitor=False) as bar:
            text_files = glob.glob(PathALL, recursive = True)
            kall = len(text_files)
        
        print(( 'Total number of '+ str(kall) + ' files were found:'+'Parsing finished! '.upper()).upper())
        ABook={}
        for i in range(len(Types)):         #Creation of Adress Book
            ABook[Types_new[i]] = []
        
# =============================================================================
#         text_files2 = [os.path.split(t)[-1] for t in text_files]
#         text_files3,Index = np.unique(text_files2,return_index=True)
#         text_files = [text_files[ii] for ii in Index]        
# ============================================================================


        for i,T in enumerate(Types_new):
            globals()['df'+ str(i)] = pd.DataFrame(ABook[T])        
        
        for i in text_files :
            
            
            for rr,Q in enumerate(i.split(os.sep)):
                temp = i.split(os.sep)[-rr].upper()
    
                Flag_anat = any([(aa in temp) for aa in T2_string])
                Flag_struct = any([(aa in temp) for aa in DTI_string])
                Flag_func = any([(aa in temp) for aa in FMRI_string])
                Flag_notAllowed = any([(aa in temp) for aa in NotAllowed])
                
                if any([Flag_anat,Flag_struct,Flag_func,Flag_notAllowed]):
                    break

            
            if not any([Flag_anat,Flag_struct,Flag_func]):
                continue
                print("The sequence type is ambiguous between registered nifti files.")
                print("To solve this problem use the following names to define sequences in ther path name:")
                print(DTI_string)
                print(FMRI_string)
                print(T2_string)
                print('Avoid using "EPI"!')
            
                
              
            if Flag_struct and not Flag_notAllowed:
                ABook["DTI"].append(i)
            if Flag_func and not Flag_notAllowed:
                ABook["rsfMRI"].append(i)
            if Flag_anat and not Flag_notAllowed:
                ABook["T2w"].append(i)

        #saving in csv file
        for n,type in enumerate(Types_new):
             if len(ABook[type]) !=0:
                 addreses= pd.DataFrame(ABook[type])
                 csv_path= "nifti_data_addreses_"+type+".csv"
                 csv_path= os.path.join(saving_path,csv_path)
                 addreses.to_csv(csv_path, sep=',',index=False)


    print('\n\ncsv files were created:' + str(saving_path))
    print('\n\n%%%%%%%%%%%%% End of the stage 1 %%%%%%%%%%%%%%%'.upper())
    print('\nStarting Stage 2 ...'.upper())
    #print('\nChosen Sequences are: ')
    #print(sequence_types)
    print('\nCalculating features...\n'.upper())
    print('This might take some time (hours/days) depending on the size of the dataset!:) ...\n\n')
    if format_type=="raw":
        fc.CheckingRawFeatures(saving_path)
        QC.toc()
    elif format_type=="nifti":
        fc.CheckingNiftiFeatures(saving_path)
        QC.toc()
    
     
    print('Plotting quality features...\n'.upper())
    QC.QCPlot(saving_path)
    QC.QCtable(saving_path)

    # remove addressed files
    for file in glob.glob(os.path.join(saving_path, '*data_addreses*.csv')) :
        os.remove(file) 
        
    # relocate calculated fearures

    calculated_features = os.path.join(saving_path, "calculated_features")
    os.mkdir(calculated_features) 
    old_files=[]
    for old_file in glob.glob(os.path.join(saving_path, '*caculated_features*.csv')) :
        old_files.append(old_file)

    new_direction= os.path.join(saving_path,"calculated_features")
    new_files=[]
    for old_file in old_files:
        path_split= os.path.split(old_file)
        files_name= path_split[1]
        join_path= os.path.join(new_direction,files_name)
        new_files.append(join_path)
                         
    for new_file,old_file in enumerate(old_files) :

        shutil.move(old_file , new_files[new_file])
    
    print('\n\n%%%%%%%%%%%%%Quality feature plots were successfully created and saved%%%%%%%%%%%%%%%\n\n'.upper())
    
    print('------------------------------------------------------------')
    print('Thank you for using our code. For questions please contact us via:')
    print('aref.kalantari-sarcheshmeh@uk-koeln.de or markus.aswendt@uk-koeln.de')
    print('Lab: AG Neuroimaging and neuroengineering of experimental stroke University Hospital Cologne')
    print('Web:https://neurologie.uk-koeln.de/forschung/ag-neuroimaging-neuroengineering/')
    print('------------------------------------------------------------')





    
    

    
    
    
    
    
    
    
    
    