import numpy as np
import os
import pandas as pd
import glob
import matplotlib.backends.backend_pdf
import nibabel as nib
import nibabel as nii
import matplotlib.pyplot as plt
import pv_conv2Nifti as pr
import alive_progress as ap
from QC import *
#%% Feature calculation of the pipeline. Core Unit of the Pipeline     
def CheckingRawFeatures(Path):   
    #Path=r"C:\Users\Erfan\Downloads\Compressed\proc_data\P5"  
    Abook = []
    Names =[]
    for file in glob.glob(os.path.join(Path, '*addreses*.csv')) :
        
        if "T2w" in file :    
             t2w_path= file
             t2w_addreses= pd.read_csv(t2w_path)
             Abook.append(t2w_addreses)
             Names.append("T2w")

        elif "DTI" in file:
            dti_path= file
            dti_addreses= pd.read_csv(dti_path)
            Abook.append(dti_addreses)
            Names.append("DTI")
        elif "fMRI" in file :
            fmri_path= file
            fmri_addreses= pd.read_csv(fmri_path)
            Abook.append(fmri_addreses)
            Names.append("rsfMRI")
      
    
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
            if kk > 0:
                print(str(kk) + 'faulty files were found:All faulty files are available in the Errorlist tab in the Excel outputs\n')
            
            print(N+' processing... \n')
            
            text_files0 = Abook2[ii]
            text_files = [i[0] for i in text_files0]
            
            dd = 1
            snrCh_vec =[]
            tsnr_vec = []
            snr_normal_vec = []
            SpatRes_vec = []
            MI_vec_all = []
            LMV_all = []
            GMV_all = []
            Max_mov_between_all = []
            text_files_new = []
            GMetric_vec =  []
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
                            input_file = nii.squeeze_image(pv.nim)
                        except ValueError:
                            ErorrList.append(tf)
                            print("Value Error: catched")
                            continue
                        except SystemError:
                            ErorrList.append(tf)
                            print("System Error: catched")
                            continue
                        
                    else:
                        ErorrList.append(tf)
                        kk = kk+1
                        continue
                   
                   
                    ########### Slice extraction 
                    selected_img = Image_Selection(input_file)                    
                    qc_path = os.path.join(saving_path,"manual_slice_inspection")
                    if not os.path.isdir(qc_path):
                        os.mkdir(qc_path)
                    img_name = str.split(tf,os.sep)[-2]
                    
                    #plt.figure()          
                    plt.axis('off')
                    plt.imshow(selected_img,cmap='gray')
                    svg_path = os.path.join(qc_path,img_name+"_"+str(N)+str(dd)+".tiff").replace(".nii","").replace(".gz","")
                    dd = dd +1
                    plt.savefig(svg_path)
                    ########### Slice extraction               
                    # other Features
                    SpatRes = ResCalculator(input_file)
                    GMetric = GoastCheck(input_file)
                    
                    if N == 'T2w':
                        
                        # Signal 2 noise ratio
                        snrCh = snrCalclualtor_chang(input_file)
                        snr_normal = snrCalclualtor_normal(input_file)   
                        
                        LMV_all = np.nan
                        GMV_all = np.nan
                        Max_mov_between_all = np.nan
                        snr_normal_vec.append(snr_normal)
                        snrCh_vec.append(snrCh)
                        
                    if N == 'DTI':
                        # Signal 2 noise ratio
                        #print(tf)
                        snrCh = snrCalclualtor_chang(input_file)
                        snr_normal = snrCalclualtor_normal(input_file)   
                        Final,Max_mov_between,GMV,LMV = Ismovement(input_file)
                        
                        
                        GMV_all.append(GMV)
                        LMV_all.append(LMV)
                        MI_vec_all.append(Final)
                        Max_mov_between_all.append(Max_mov_between)
                        snr_normal_vec.append(snr_normal)
                        snrCh_vec.append(snrCh)
                        
                    if N == 'rsfMRI':
                        #temporal signal 2 noise ratio
                        #print(tf)
                        tSNR = TsnrCalclualtor(input_file)
                        Final,Max_mov_between,GMV,LMV = Ismovement(input_file)
                        
                        Max_mov_between_all.append(Max_mov_between)
                        GMV_all.append(GMV)
                        LMV_all.append(LMV)
                        MI_vec_all.append(Final)
                        Max_mov_between_all.append(Max_mov_between)
                        tsnr_vec.append(tSNR)

                        
                        
                    
                    i=i+1
                    bar()
                    text_files_new.append(tf)
                    SpatRes_vec.append(SpatRes)
                    GMetric_vec.append(GMetric)
                    
            # Saving parsed files to excel sheets
            #AR = [text_files_new,np.array(SpatRes_vec),np.array(snrCh_vec),np.array(LMV_all),np.array(GMV_all),np.array(snr_normal_vec)]        
            
            # using the savetxt 
            # from the numpy module
            
            df = pd.DataFrame()
            df['FileAddress'] = text_files_new
            df['SpatRx'] = np.array(SpatRes_vec)[:,0]
            df['SpatRy'] = np.array(SpatRes_vec)[:,1]
            df['Slicethick'] = np.array(SpatRes_vec)[:,2]
            df['Goasting'] = np.array(GMetric_vec)
            
            
            if N == 'T2w':
                 df['SNR Chang'] = np.array(snrCh_vec)
                 df['SNR Normal'] = np.array(snr_normal_vec)
                 
            elif N == 'DTI':
                 df['SNR Chang'] = np.array(snrCh_vec)
                 df['SNR Normal'] = np.array(snr_normal_vec)
                 df['Displacement factor (std of Mutual information)']=np.array(LMV_all)
                 #df['Maximal displacement']=AR[4]
                 
            elif N == "rsfMRI":
                 df['tSNR (Averaged Brain ROI)'] = np.array(tsnr_vec)
                 df['Displacement factor (std of Mutual information)']=np.array(LMV_all)
                 #df['Maximal displacement']=AR[4]
                 
            if N=="T2w":
                t2w_result= os.path.join(Path,"caculated_features_T2w.csv")
                df.to_csv( t2w_result)

            elif N=="DTI":    
                dti_result= os.path.join(Path,"caculated_features_DTI.csv")
                df.to_csv( dti_result)   

            elif N=="rsfMRI":
                fmri_result= os.path.join(Path,"caculated_features_fMRI.csv")
                df.to_csv(fmri_result)
        
          
    
    print('\n\noutput files was created:' + str(Path))
    
    print('\n\n%%%%%%%%%%%%%End of the Second stage%%%%%%%%%%%%%%%\n\n'.upper())
   
#%% exact above function but this time for nifti format
def CheckingNiftiFeatures(Path):   
    
    
    Abook = []
    Names =[]
    for file in glob.glob(os.path.join(Path, '*addreses*.csv')) :

        if "T2w" in file :    
             t2w_path= file
             t2w_addreses= pd.read_csv(t2w_path)
             Abook.append(t2w_addreses)
             Names.append("T2w")             
        elif "DTI" in file:
            dti_path= file
            dti_addreses= pd.read_csv(dti_path)
            Abook.append(dti_addreses)
            Names.append("DTI")
        elif "fMRI" in file :
            fmri_path= file
            fmri_addreses= pd.read_csv(fmri_path)
            Abook.append(fmri_addreses)
            Names.append("rsfMRI")
    
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
            
            dd = 1
            snrCh_vec =[]
            tsnr_vec = []
            snr_normal_vec = []
            SpatRes_vec = []
            MI_vec_all = []
            LMV_all = []
            GMV_all = []
            Max_mov_between_all = []
            text_files_new = []
            GMetric_vec =  []
            kk = 0
            i=1
            
            
            with ap.alive_bar(len(text_files),spinner='wait') as bar:

                for tf in text_files:

                    tf = str(tf)
                   
                    if "DTI".upper() in tf.upper() :
                        N= "DTI"
                    if  "T2w".upper() in tf.upper():
                        N="T2w"
                    if  "fMRI".upper() in tf.upper():
                        N="rsfMRI"
                    tf = os.path.normpath(tf)
                    input_file= nib.load(tf)
                   
                    ########### Slice extraction 
                    selected_img = Image_Selection(input_file)                    
                    qc_path = os.path.join(Path,"manual_slice_inspection")
                    if not os.path.isdir(qc_path):
                        os.mkdir(qc_path)
                    img_name = str.split(tf,os.sep)[-1]
                    #plt.figure()          
                    plt.axis('off')
                    plt.imshow(selected_img,cmap='gray')
                    svg_path = os.path.join(qc_path,img_name+"_"+str(N)+str(dd)+".tiff").replace(".nii","").replace(".gz","")
                    dd = dd +1
                    plt.savefig(svg_path)
                    ########### Slice extraction               
                    # other Features
                    SpatRes = ResCalculator(input_file)
                    GMetric = GoastCheck(input_file)
                    
                    if N == 'T2w':
                        # Signal 2 noise ratio
                        snrCh = snrCalclualtor_chang(input_file)
                        snr_normal = snrCalclualtor_normal(input_file)   
                        
                        LMV_all = np.nan
                        GMV_all = np.nan
                        Max_mov_between_all = np.nan
                        snr_normal_vec.append(snr_normal)
                        snrCh_vec.append(snrCh)
                        
                    if N == 'DTI':
                        # Signal 2 noise ratio
                        
                        snrCh = snrCalclualtor_chang(input_file)
                        snr_normal = snrCalclualtor_normal(input_file)   
                        Final,Max_mov_between,GMV,LMV = Ismovement(input_file)
                        
                        
                        GMV_all.append(GMV)
                        LMV_all.append(LMV)
                        MI_vec_all.append(Final)
                        Max_mov_between_all.append(Max_mov_between)
                        snr_normal_vec.append(snr_normal)
                        snrCh_vec.append(snrCh)
                        
                    if N == 'rsfMRI':
                        #temporal signal 2 noise ratio
                        tSNR = TsnrCalclualtor(input_file)
                        Final,Max_mov_between,GMV,LMV = Ismovement(input_file)
                        
                        Max_mov_between_all.append(Max_mov_between)
                        GMV_all.append(GMV)
                        LMV_all.append(LMV)
                        MI_vec_all.append(Final)
                        Max_mov_between_all.append(Max_mov_between)
                        tsnr_vec.append(tSNR)
                    
                    
                    
                    i=i+1
                    bar()
                    text_files_new.append(tf)
                    SpatRes_vec.append(SpatRes)
                    GMetric_vec.append(GMetric)
                    
                     
                     
            # Saving parsed files to excel sheets
            #AR = [text_files_new,np.array(SpatRes_vec),np.array(snrCh_vec),np.array(LMV_all),np.array(GMV_all),np.array(snr_normal_vec)]
            
            
            # using the savetxt 
            # from the numpy module
            
            df = pd.DataFrame()
            df['FileAddress'] = text_files_new
            df['SpatRx'] = np.array(SpatRes_vec)[:,0]
            df['SpatRy'] = np.array(SpatRes_vec)[:,1]
            df['Slicethick'] = np.array(SpatRes_vec)[:,2]
            df['Goasting'] = np.array(GMetric_vec)
            
          
            
            if N == 'T2w':
                 df['SNR Chang'] = np.array(snrCh_vec)
                 df['SNR Normal'] = np.array(snr_normal_vec)
                 
            elif N == 'DTI':
                 df['SNR Chang'] = np.array(snrCh_vec)
                 df['SNR Normal'] = np.array(snr_normal_vec)
                 df['Displacement factor (std of Mutual information)']=np.array(LMV_all)
                 #df['Maximal displacement']=AR[4]
                 
            elif N == "rsfMRI":
                 df['tSNR (Averaged Brain ROI)'] = np.array(tsnr_vec)
                 df['Displacement factor (std of Mutual information)']=np.array(LMV_all)
                 #df['Maximal displacement']=AR[4]
                 
            if N=="T2w":
                t2w_result= os.path.join(Path,"caculated_features_T2w.csv")
                df.to_csv( t2w_result)

            elif N=="DTI":    
                dti_result= os.path.join(Path,"caculated_features_DTI.csv")
                df.to_csv( dti_result)   

            elif N=="rsfMRI":
                fmri_result= os.path.join(Path,"caculated_features_fMRI.csv")
                df.to_csv(fmri_result)
        
    
    print('\n\noutput file was created:' + str(Path))
    
    print('\n\n%%%%%%%%%%%%%End of the Second stage%%%%%%%%%%%%%%%\n\n'.upper())
    
   
#%%