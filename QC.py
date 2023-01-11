#%% 
"""
Version 1.0
Name: Aref Kalantari
Email: aref.kalantari-sarcheshmeh@uk-koeln.de
Date: 24.08.2021 - 02.03.2022
-----------------------------
Code Describtion: Quality Control Toolbox. Every tool (function) needed can be found here and be modified.
-----------------------------
Lab: AG Neuroimaging and neuroengineering of experimental stroke 
Supervisor: Dr. rer. nat. Markus Aswendt (markus.aswendt@uk-koeln.de)
"""

#%% Loading nececcery libraries

import numpy as np
import changSNR as ch
from matplotlib.pyplot import imshow
import os
import pandas as pd
import glob
from openpyxl import load_workbook
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
import seaborn as sns
import matplotlib.patches as mpatches
import numpy as np
import openpyxl
import nibabel as nib
import glob
import pv_conv2Nifti as pr
import openpyxl
import nibabel as nii
import alive_progress as ap
import time
#%% Res function

def ResCalculator(input_file):
    
    HDR = input_file.header
    Spati_Res = HDR['pixdim'][1:4]
    
    return Spati_Res



#%% SNR function
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
    #for slc in range(ns_lower,ns_upper):
    for slc in range(ns_lower,ns_upper):    
        #   Print % of progress
        #print('S' + str(slc + 1), end=",")

        # Decision if the input data is DTI type or T2w
        if nd == 3:
            slice = imgData[:, :, slc]
           
        if nd == 4:
            slice = imgData[:, :, slc,10]
            
        curSnrCHMap, estStdChang, estStdChangNorm = ch.calcSNR(slice, 0, 1)
        noiseChSNR[slc] = estStdChang

    snrCh = 20 * np.log10(np.mean(imgData) / np.mean(noiseChSNR))

    return snrCh

#%% SNR function
def snrCalclualtor_nifti(input_file):

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
            slice = imgData[:, :, slc]
            
        curSnrCHMap, estStdChang, estStdChangNorm = ch.calcSNR(slice, 0, 1)
        noiseChSNR[slc] = estStdChang

    snrCh = 20 * np.log10(np.mean(imgData) / np.mean(noiseChSNR))

    return snrCh


#%% TSNR function

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


#%% Calculating Mutual Information: based on https://matthew-brett.github.io/teaching/mutual_information.html


def mutualInfo(Im1,Im2):

    t1_slice = Im1
    t2_slice = Im2

    hist_2d, x_edges, y_edges = np.histogram2d(t1_slice.ravel(),t2_slice.ravel(),bins=20)

    hist_2d_log = np.zeros(hist_2d.shape)
    non_zeros = hist_2d != 0
    hist_2d_log[non_zeros] = np.log(hist_2d[non_zeros])
    
    pxy = hist_2d / float(np.sum(hist_2d))
    px = np.sum(pxy, axis=1) # marginal for x over y
    py = np.sum(pxy, axis=0) # marginal for y over x
    px_py = px[:, None] * py[None, :] # Broadcast to multiply marginals
    # Now we can do the calculation using the pxy, px_py 2D arrays
    nzs = pxy > 0 # Only non-zero pxy values contribute to the sum
    MI = np.sum(pxy[nzs] * np.log(pxy[nzs] / px_py[nzs]))
    
    return MI


#%% Movement detection of rsFRI function (based on mutual information)

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
        
        MI = mutualInfo(Im_fix,Im_rot[:,:,z])
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

#%% Getting range 

def getrange(numbers):
    return max(numbers) - min(numbers)


#%% Plotting QC Histogram and etc.

def QCPlot(Path):
    Names = []
    xls = pd.ExcelFile(Path,engine= 'openpyxl')
    Names = xls.sheet_names
    saving_path = os.path.dirname(Path) 
    
    if 'ErrorData' in Names:
        Names.remove('ErrorData')
    
    Abook = []
    for n in Names:
        Abook.append(pd.read_excel(Path,engine= 'openpyxl',sheet_name = n))
    
    # creating plots
    sns.set_palette("colorblind")
    #f1, ax1 = plt.subplots(3,3)
    rr = 1
    hh = 1
    
    
    for nn,N in enumerate(Names):
        COL = list(Abook[nn].columns)
        COL.pop(0)
        D = Abook[nn]
        for cc,C in enumerate(COL):
            Data = list(D[C])
            if C== 'SNR Chang' or C == 'tSNR Chang' or C=='Local Movement Variability':
                #plot histogrm
                plt.figure(hh,figsize=(10,5))
                ax2 = plt.subplot(1,1,1)
                for dd,DD in enumerate(Data):
                    if DD == np.inf:
                        Data[dd] = np.nan
                
                q75, q25 = np.nanpercentile(Data, [75 ,25])
                iqr = q75 - q25
                
                B = round((np.nanmax(Data)-np.nanmin(Data)) / (2 * iqr / (len(Data)**(1/3))))
               
                y, x, bars = plt.hist(Data, bins= B, histtype= 'bar',edgecolor='white')
                plt.xlabel(N+': '+C + ' [a.u.]')
                plt.ylabel("Frequency")
                ax2.spines['right'].set_visible(False)
                ax2.spines['top'].set_visible(False)
                #calculate interquartile range of values in the 'points' column
                
                if C == 'Local Movement Variability':
                    ll = q75+1.5*iqr
                    plt.text(1.07*ll,2*max(y)/3,'Q3 + 1.5*IQ',color='grey')
                    for b,bar in enumerate(bars):
                        if bar.get_x() > ll:
                          bar.set_facecolor("red")
                          
                    
                else:
                    ll = q25-1.5*iqr
                    plt.text(1.001*ll,2*max(y)/3,'Q1 - 1.5*IQ',color='grey')
                    for b,bar in enumerate(bars):
                        if bar.get_x() < ll:
                            bar.set_facecolor("red")
                            
                plt.axvline(ll,color = 'grey',linestyle='--')
                plt.suptitle(N+': '+C,weight="bold")
                hh = hh + 1
                
    
                red_patch = mpatches.Patch(color='red', label='Discard')
                blue_patch = mpatches.Patch(color='tab:blue', label='Keep')
                plt.legend(handles=[blue_patch,red_patch])
                   # plt.savefig(os.path.dirname(Path) + "\ResHomogenity.png",dpi=300)
               
               
    plt.figure(hh,figsize=(14, 10))
    for nn,N in enumerate(Names):
        COL = list(Abook[nn].columns)
        COL.pop(0)
        D = Abook[nn]
        for cc,C in enumerate(COL):
            Data = list(D[C])           
            if C == 'SpatRx' or C == 'SpatRy' or C == 'Slicethick':
                #plot pieplots
                labels = list(set(Data))
                sizes = [Data.count(l) for l in labels]
                labels= list(np.round(labels,3))
                labels2=[str(l)+' mm' for l in labels]
                
                ax1 = plt.subplot(len(Names),3,rr)           
                ax1.pie(sizes, labels=labels2, autopct='%1.0f%%', startangle=90)
                ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
                ax1.set_title(N+':'+C)
                plt.suptitle('Resolution homogeneity between data',weight="bold")
                rr = rr+1
    
    pdf = matplotlib.backends.backend_pdf.PdfPages(os.path.dirname(Path) + "/QCfigures.pdf")
    for fig in range(1, plt.figure().number): ## will open an empty extra figure :(
        pdf.savefig( fig )
    pdf.close()

#%% Adjusting the existing feature table by adding a new sheet to it with the data that need to be discarded

def QCtable(Path):
    

    Abook = []
    Names =[]
    for file in glob.glob(Path+ '*caculated_features*.csv') :
        
        if "DTI" in file:
            dti_path= file
            dti_features= pd.read_csv(dti_path)
            Abook.append(dti_features)
            Names.append("DTI")
        elif "fMRI" in file :
            fmri_path= file
            fmri_features= pd.read_csv(fmri_path)
            Abook.append(fmri_features)
            Names.append("rsfMRI")
        elif "T2w" in file :    
             t2w_path= file
             t2w_features= pd.read_csv(t2w_path)
             Abook.append(t2w_features)
             Names.append("T2w")    

    
    
  
    ST = []
    COE = []
    AvV = []
    V = []
    Pathes = []
    Med = []
    MaX = []
    MiN= []
    for nn,N in enumerate(Names):

            
        d= Abook[nn]
        COL = Abook[nn].columns
        
        for cc,C in enumerate(COL):
            
            D = d[C]
            
            
            if C == 'SNR Chang' or C == 'tSNR Chang':
                
                for dd,DD in enumerate(D):
                    
                    if DD == np.inf:
                        D[dd] = np.nan
                        
                
                q75, q25 = np.nanpercentile(D, [75 ,25])
                
                iqr = q75 - q25
                ll = q25-1.5*iqr #lower limit
                Index = D<ll
                
                P = d[COL[1]][Index]
                
                M = D.mean()
                Me = D.median()
                Mi = D.min()
                Ma = D.max()
              
         
                Pathes.extend(P)
                ST.extend([N]*len(P))
                COE.extend([C]*len(P))
                AvV.extend([M]*len(P))
                V.extend(D[Index])
                Med.extend([Me]*len(P))
                MiN.extend([Mi]*len(P))
                MaX.extend([Ma]*len(P))
                
                
            if C == 'Local Movement Variability':
                q75, q25 = np.nanpercentile(D, [75 ,25])
                iqr = q75 - q25
                ul = q75+1.5*iqr #upper limit
                Index = D>ul
                P = d[COL[1]][Index]
                M = D.mean()
                Me = D.median()
                Mi = D.min()
                Ma = D.max()
                 
                Pathes.extend(P)
                ST.extend([N]*len(P))
                COE.extend([C]*len(P))
                AvV.extend([M]*len(P))
                V.extend(D[Index])
                Med.extend([Me]*len(P))
                MiN.extend([Mi]*len(P))
                MaX.extend([Ma]*len(P))
                
    
            if N == 'ErrorData': 
                Pathes.extend(D)
                S = 'Faulty Data'
                ST.extend([S]*len(D))
                COE.extend(['-']*len(D))
                AvV.extend(['-']*len(D))
                V.extend('-'*len(D))
                Med.extend('-'*len(D))
                MiN.extend('-'*len(D))
                MaX.extend('-'*len(D))
                        
    List = {"Pathes":Pathes,"Sequence Type":ST, "Problematic Quality Feature":COE,"Value":V,"Mean":AvV,"Median":Med,"Min":MiN,"Max":MaX}
    
    df = pd.DataFrame(List)
    final_result = os.path.join(Path,"unaccountable_data.csv")
    df.to_csv( final_result)    
    

    
    

#%% Feature calculation of the pipeline. Core Unit of the Pipeline
     
def CheckingrawFeatures(Path):   
    #Path=r"C:\Users\Erfan\Downloads\Compressed\proc_data\P5"  
    Abook = []
    Names =[]
    for file in glob.glob(os.path.join(Path, '*addreses*.csv')) :
        
        if "DTI" in file:
            dti_path= file
            dti_addreses= pd.read_csv(dti_path)
            Abook.append(dti_addreses)
            Names.append("DTI")
        elif "fMRI" in file :
            fmri_path= file
            fmri_addreses= pd.read_csv(fmri_path)
            Abook.append(fmri_addreses)
            Names.append("rsfMRI")
        elif "T2w" in file :    
             t2w_path= file
             t2w_addreses= pd.read_csv(t2w_path)
             Abook.append(t2w_addreses)
             Names.append("T2w")
      
    
    ErorrList = []
 
    
    saving_path = Path
    

    C = np.array([not Check.empty for Check in Abook])
    Names = np.array(Names)[C].tolist()
    Names.append('ErrorData')
    Abook2 = [a.values.tolist() for a in Abook]
    #% Calculate all the SNR for all the data that was found 
    # in the last step and saving it into a vector
    # Load Bruker data from each address 

    kk =0
    for ii,N in enumerate(Names):
        if N != 'ErrorData':
            
            print(N+' processing... \n')
            
            text_files0 = Abook2[ii]
            text_files = [i[0] for i in text_files0]
            
            
            snrCh_vec =[]
            SpatRes_vec = []
            MI_vec_all = []
            LMV_all = []
            TypeMov_all = []
            text_files_new = []
            kk = 0
            i=1
           
            with ap.alive_bar(len(text_files),spinner='wait') as bar:
                for tf in text_files:
                    
                    tf = str(tf)
                    tf = os.path.normpath(tf)
                   
                    path_split = tf.split(os.sep)
                    
                    procno = str(1)
                    expno =path_split[-1]
                    
                    study = path_split[-2]
                    raw_folder = os.sep.join(path_split[:-2])
                    proc_folder = os.path.join(raw_folder,'proc_data') #Here still adjustment is needed
                    pv = pr.Bruker2Nifti(study, expno, procno, raw_folder, proc_folder, ftype='NIFTI_GZ')
                    CP_v = os.path.join(tf ,'pdata','1','visu_pars') # Check Parameter: Visu_pars
                    CP_a = os.path.join(tf , 'acqp') # Check Parameter: acqp
                    
                    if os.path.isfile(CP_v) and os.path.isfile(CP_a):
                        try:
                            pv.read_2dseq(map_raw=False, pv6=False)
                        except ValueError:
                            ErorrList.append(tf)
                            print("Value Error: catched")
                            continue
                        except SystemError:
                            ErorrList.append(tf)
                            print("System Error: catched")
                            continue
                        
                        input_file = nii.squeeze_image(pv.nim)
                        
                    else:
                        ErorrList.append(tf)
                        kk = kk+1
                        continue
                   
                    # Resoultution
                    
                    SpatRes = ResCalculator(input_file)
                    
                    
                    if N == 'T2w' or N == 'DTI':
                        # Signal 2 noise ratio
                        try:
                            snrCh = snrCalclualtor(input_file)
                        except Exception:
                            ErorrList.append(tf)
                            kk = kk+1
                            continue
                            
                        LMV_all = np.nan
                        TypeMov_all = np.nan
                        
                    if N == 'rsfMRI':
                        #temporal signal 2 noise ratio
                        try:
                            snrCh = TsnrCalclualtor(input_file)
                        except Exception:
                            ErorrList.append(tf)
                            kk = kk+1
                            continue
                        
                        # movement severity with the help of mutual information
                        Final,GMV,LMV,TypeMov = Ismovement(input_file)
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
                 
            if N=="T2w":

            
                 
                t2w_result= os.path.join(Path,"caculated_features_T2w.csv")
                df.to_csv( t2w_result)

            elif N=="DTI":    
                dti_result= os.path.join(Path,"caculated_features_DTI.csv")
                df.to_csv( dti_result)   



            elif N=="fMRI":
                fmri_result= os.path.join(Path,"caculated_features_fMRI.csv")
                df.to_csv( fmri_result)                              
            else:
                df = pd.DataFrame()
                
                df['ErorrList'] = ErorrList
                ErorrList_result= os.path.join(Path,"Faulty_Datasets.csv")
                df.to_csv(ErorrList_result)
        
    
    print('\n\noutput files was created:' + str(Path))
    
    print('\n\n%%%%%%%%%%%%%End of the Second stage%%%%%%%%%%%%%%%\n\n'.upper())
    print('Plotting quality features...\n'.upper())
    
    #QCPlot(saving_path2)
    QCtable(Path)
    print('\n\n%%%%%%%%%%%%%Quality feature plots were successfully created and saved%%%%%%%%%%%%%%%\n\n'.upper())




#exact above function but this time for nifti format
def CheckingNiftiFeatures(Path):   
    
    
    Abook = []
    Names =[]
    for file in glob.glob(os.path.join(Path, '*addreses*.csv')) :
       
        if "DTI" in file:
            dti_path= file
            dti_addreses= pd.read_csv(dti_path)
            Abook.append(dti_addreses)
            Names.append("DTI")
        elif "fMRI" in file :
            fmri_path= file
            fmri_addreses= pd.read_csv(fmri_path)
            Abook.append(fmri_addreses)
            Names.append("rsfMRI")
        elif "T2w" in file :    
             t2w_path= file
             t2w_addreses= pd.read_csv(t2w_path)
             Abook.append(t2w_addreses)
             Names.append("T2w")
    
    ErorrList = []
 
    
    saving_path = os.path.dirname(Path) 
    

    C = np.array([not Check.empty for Check in Abook])
    Names = np.array(Names)[C].tolist()
    Names.append('ErrorData')
    Abook2 = [a.values.tolist() for a in Abook]


    #% Calculate all the SNR for all the data that was found 
    # in the last step and saving it into a vector
    # Load Bruker data from each address 

    kk =0
    
    for ii,N in enumerate(Names):
        
        if N != 'ErrorData':
            if kk > 0:
                print(str(kk) + 'faulty files were found:All faulty files are available in the Errorlist tab in the Excel outputs\n')
            
            print(N+' processing... \n')
            
            text_files0 = Abook2[ii]
            text_files = [i[0] for i in text_files0]
            
            
            snrCh_vec =[]
            SpatRes_vec = []
            MI_vec_all = []
            LMV_all = []
            TypeMov_all = []
            text_files_new = []
            kk = 0
            i=1
            
            
            with ap.alive_bar(len(text_files),spinner='wait') as bar:

                for tf in text_files:

                    tf = str(tf)
                   
                    if "DTI" in tf :
                        N= "DTI"
                    if  "T2w" in tf:
                        N="T2w"
                    if  "fMRI" in tf:
                        N="fMRI"
                    tf = os.path.normpath(tf)
                    
                    

                    input_file= nib.load(tf)
                   
                    # Resoultution
                    SpatRes = ResCalculator(input_file)
                    
                    
                    
                    
                    if N == 'T2w' or N == 'DTI':
                        # Signal 2 noise ratio
                        
                        
                        snrCh = snrCalclualtor_nifti(input_file)
                           
                        
                            
                            
                        
                            
                        LMV_all = np.nan
                        TypeMov_all = np.nan
                        
                    if N == 'fMRI':
                        #temporal signal 2 noise ratio
                        
                        snrCh = TsnrCalclualtor(input_file)
                        
                            
                        
                        # movement severity with the help of mutual information
                        Final,GMV,LMV,TypeMov = Ismovement(input_file)
                        TypeMov_all.append(TypeMov)
                        LMV_all.append(LMV)
                        MI_vec_all.append(Final)


                    snrCh_vec.append(snrCh)
                    
                    i=i+1
                    text_files_new.append(tf)
                    bar()
                    SpatRes_vec.append(SpatRes) 
                     
                
                     
                     
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
                 
            if N=="T2w":

            
                
                t2w_result= os.path.join(Path,"caculated_features_T2w.csv")
                df.to_csv( t2w_result)

            elif N=="DTI":    
                dti_result= os.path.join(Path,"caculated_features_DTI.csv")
                df.to_csv( dti_result)   



            elif N=="fMRI":
                fmri_result= os.path.join(Path,"caculated_features_fMRI.csv")
                df.to_csv( fmri_result)                              
            else:
                df = pd.DataFrame()
                
                df['ErorrList'] = ErorrList
                df.to_csv( r"C:\Users\Erfan\Downloads\Compressed\proc_data\P5\Sham\df.csv")
        
    
    print('\n\nExcel file was created:' + str(Path))
    
    print('\n\n%%%%%%%%%%%%%End of the Second stage%%%%%%%%%%%%%%%\n\n'.upper())
    print('Plotting quality features...\n'.upper())
    QCtable(Path)
    #QCPlot(dti_result,fmri_result,t2w_result)
    
    print('\n\n%%%%%%%%%%%%%Quality feature plots were successfully created and saved%%%%%%%%%%%%%%%\n\n'.upper())
#%% Tic Toc Timer


def TicTocGenerator():
    # Generator that returns time differences
    ti = 0           # initial time
    tf = time.time() # final time
    while True:
        ti = tf
        tf = time.time()
        yield tf-ti # returns the time difference

TicToc = TicTocGenerator() # create an instance of the TicTocGen generator

# This will be the main function through which we define both tic() and toc()
def toc(tempBool=True):
    # Prints the time difference yielded by generator instance TicToc
    tempTimeInterval = next(TicToc)
    if tempBool:
        print( "Elapsed time: %f seconds.\n" %tempTimeInterval )

def tic():
    # Records a time in TicToc, marks the beginning of a time interval
    toc(False)
#%% For Questions please Contact: aref.kalantari-sarcheshmeh@uk-koeln.de