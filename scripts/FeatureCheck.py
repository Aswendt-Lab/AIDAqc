import numpy as np
import os
import pandas as pd
import glob
import matplotlib.backends.backend_pdf
import nibabel as nib
from nilearn.plotting import plot_img
from nilearn.image import index_img
import pv_conv2Nifti as pr
import alive_progress as ap
import pv_parser as par
from QC import *
#%% Feature calculation of the pipeline. Core Unit of the Pipeline     
def CheckingRawFeatures(Path):
    #Path=r"C:\Users\Erfan\Downloads\Compressed\proc_data\P5"  
    Abook = []
    Names =[]
    for file in glob.glob(os.path.join(Path, '*addreses*.csv')) :
        
        if "anat" in file :    
             t2w_path= file
             t2w_addreses= pd.read_csv(t2w_path)
             Abook.append(t2w_addreses)
             Names.append("anat")

        elif "diff" in file:
            dti_path= file
            dti_addreses= pd.read_csv(dti_path)
            Abook.append(dti_addreses)
            Names.append("diff")
        elif "func" in file :
            fmri_path= file
            fmri_addreses= pd.read_csv(fmri_path)
            Abook.append(fmri_addreses)
            Names.append("func")
    
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
                print(str(kk) + ' 1 faulty file(s) were found. Note: these files will be listed in CanNotProcessTheseFiles.csv \n')
            
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
            img_names_new = []
            keys = []
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
                            input_file = nib.squeeze_image(pv.nim)
                            # this is a workaround to convert the data for nilearn plotting
                            affine = np.eye(4)
                            affine[0, 0] = pv.nim.header.get('pixdim')[1]
                            affine[1, 1] = pv.nim.header.get('pixdim')[2]
                            affine[2, 2] = pv.nim.header.get('pixdim')[3]
                            affine[3, 3] = pv.nim.header.get('pixdim')[4]
                            input_file=nib.Nifti1Image(input_file.get_fdata(), affine=affine,dtype=int(32))
                        except ValueError:
                            ErorrList.append(tf+"_Value Error")
                            print(tf)
                            print("Value Error: catched")
                            continue
                        except SystemError:
                            ErorrList.append(tf+"_System Error")
                            print(tf)
                            print("System Error: catched")
                            continue
                        except KeyError:
                            ErorrList.append(tf+"_KeyError")
                            print(tf)
                            print("KeyError: catched")
                            continue
                        except FileNotFoundError: 
                            ErorrList.append(tf+"_FileNotFoundError")
                            print(tf)
                            print("System Error: catched")
                            continue
                        
                    else:
                        ErorrList.append(tf+"No_Visu_pars_file")
                        print("No Visu_pars file found")
                        kk = kk+1
                        continue

                    # determine sequence name 
                    NameTemp = par.read_param_file(CP_a)
                    MN = NameTemp[1]["ACQ_method"].upper()  #Here we check what the name of the sequence is
                    MN2 = NameTemp[1]["ACQ_protocol_name"].upper()
                    KEY = MN + MN2
                    keys.append(KEY)
                   
# =============================================================================
#                     Size = ((input_file.get_fdata()).nbytes/(1024*1024))
#                     if Size < 3:
#                         ErorrList.append(tf)
#                         continue
# =============================================================================
                    ########### Slice extraction 
                    qc_path = os.path.join(saving_path,"manual_slice_inspection")
                    if not os.path.isdir(qc_path):
                        os.mkdir(qc_path)
                    img_name = str.split(tf,os.sep)[-2]
                    full_img_name = str(N)+"_" + img_name+"_"+ str(dd)+".png".replace(".nii","").replace(".gz","")
                    img_names_new.append(full_img_name) 
                    
                    # grandjean patch to output ortho representation of the image in manual slice inspection
                    svg_path = os.path.join(qc_path,str(N)+"_"+ img_name+"_"+ str(dd)+".png").replace(".nii","").replace(".gz","")
                    dd = dd +1
                    if len(input_file.shape) == 4:
                        input_file_img = index_img(input_file,0)
                    else:
                        input_file_img = input_file
                    plot_img(input_file_img, title=full_img_name, output_file=svg_path)
                    ########### Slice extraction               
                    # other Features
                    SpatRes = ResCalculator(input_file)
                    GMetric = GhostCheck(input_file)
                    
                    
                    if N == 'anat':
                        
                        # Signal 2 noise ratio
                        snrCh = snrCalclualtor_chang(input_file)
                        snr_normal = snrCalclualtor_normal(input_file)   
                        
                        LMV_all = np.nan
                        GMV_all = np.nan
                        Max_mov_between_all = np.nan
                        snr_normal_vec.append(snr_normal)
                        snrCh_vec.append(snrCh)
                        
                    elif N == 'diff':
                        # Signal 2 noise ratio
                        #print(tf)
                        snrCh = snrCalclualtor_chang(input_file)
                        snr_normal = snrCalclualtor_normal(input_file)   
                        Final,Max_mov_between,GMV,LMV = Ismotion(input_file)
                        
                        
                        GMV_all.append(GMV)
                        LMV_all.append(LMV)
                        MI_vec_all.append(Final)
                        Max_mov_between_all.append(Max_mov_between)
                        snr_normal_vec.append(snr_normal)
                        snrCh_vec.append(snrCh)
                        
                    elif N == 'func':
                        #temporal signal 2 noise ratio
                        #print(tf)
                        tSNR = TsnrCalclualtor(input_file)
                        Final,Max_mov_between,GMV,LMV = Ismotion(input_file)
                        
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
            
            df = pd.DataFrame()
            df['FileAddress'] = text_files_new
            df['sequence name'] = keys
            df["corresponding_img"] = img_names_new
            df['SpatRx'] = np.array(SpatRes_vec)[:,0]
            df['SpatRy'] = np.array(SpatRes_vec)[:,1]
            df['Slicethick'] = np.array(SpatRes_vec)[:,2]
            df['Ghosting'] = np.array(GMetric_vec)
            
            
            if N == 'anat':
                 df['SNR Chang'] = np.array(snrCh_vec)
                 df['SNR Normal'] = np.array(snr_normal_vec)
                 
            elif N == 'diff':
                 df['SNR Chang'] = np.array(snrCh_vec)
                 df['SNR Normal'] = np.array(snr_normal_vec)
                 df['Displacement factor (std of Mutual information)']=np.array(LMV_all)
                 #df['Maximal displacement']=AR[4]
                 
            elif N == "func":
                 df['tSNR (Averaged Brain ROI)'] = np.array(tsnr_vec)
                 df['Displacement factor (std of Mutual information)']=np.array(LMV_all)
                 #df['Maximal displacement']=AR[4]
                 
            if N=="anat":
                t2w_result= os.path.join(Path,"caculated_features_anat.csv")
                df.to_csv( t2w_result)

            elif N=="diff":    
                dti_result= os.path.join(Path,"caculated_features_diff.csv")
                df.to_csv( dti_result)   

            elif N=="func":
                fmri_result= os.path.join(Path,"caculated_features_func.csv")
                df.to_csv(fmri_result)

    if ErorrList:            
        dfNewError = pd.DataFrame(ErorrList)
        new_file_path = os.path.join(saving_path, "CanNotProcessTheseFiles.csv")
        dfNewError.to_csv(new_file_path)
      
    print('\n\noutput files were created:' + str(Path))
    
    print('\n\n%%%%%%%%%%%%% End of the stage 2 %%%%%%%%%%%%%%%\n\n'.upper())
   
#%% exact above function but this time for nifti format
def CheckingNiftiFeatures(Path):   
    
    
    Abook = []
    Names =[]
    for file in glob.glob(os.path.join(Path, '*addreses*.csv')) :

        if "anat" in file :    
             t2w_path= file
             t2w_addreses= pd.read_csv(t2w_path)
             Abook.append(t2w_addreses)
             Names.append("anat")             
        elif "diff" in file:
            dti_path= file
            dti_addreses= pd.read_csv(dti_path)
            Abook.append(dti_addreses)
            Names.append("diff")
        elif "func" in file :
            fmri_path= file
            fmri_addreses= pd.read_csv(fmri_path)
            Abook.append(fmri_addreses)
            Names.append("func")
    
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
                print(str(kk) + 'faulty files were found: All faulty files are available in the Errorlist tab in the Excel outputs\n')
            
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
            img_names_new = []
            GMetric_vec =  []
            kk = 0
            i=1
            
            
            with ap.alive_bar(len(text_files),spinner='wait') as bar:

                for tf in text_files:

                    tf = str(tf)
                    print(tf)
                    """
                    if "DTI".upper() in tf.upper() :
                            N= "DTI"
                    if  "T2w".upper() in tf.upper():
                            N="T2w"
                    if  "fMRI".upper() in tf.upper():
                            N="rsfMRI"
                    """
                    tf = os.path.normpath(tf)
                    try:
                        input_file= nib.load(tf)
                    except nib.loadsave.ImageFileError:
                        print("could not load the following file (check the size of the file):")
                        print(tf)
                        ErorrList.append(tf)
                        continue
                    
                    ########### Slice extraction 
                    qc_path = os.path.join(Path,"manual_slice_inspection")
                    if not os.path.isdir(qc_path):
                        os.mkdir(qc_path)
                    img_name = str.split(tf,os.sep)[-1]
                    folder_name = str.split(tf,os.sep)[-2] 
                    full_img_name = (str(N)+"_"+folder_name+"_"+img_name+"_"+str(dd)+".png").replace(".nii","").replace(".gz","")
                    img_names_new.append(full_img_name)
                    
                    # grandjean patch to output ortho representation of the image in manual slice inspection
                    svg_path = os.path.join(qc_path,str(N)+"_"+ img_name+"_"+ str(dd)+".png").replace(".nii","").replace(".gz","")
                    dd = dd +1
                    if len(input_file.shape) == 4:
                        input_file_img = index_img(input_file,0)
                    else:
                        input_file_img = input_file
                    plot_img(input_file_img, title=full_img_name, output_file=svg_path)
                    ########### Slice extraction               
                    # other Features
                    SpatRes = ResCalculator(input_file)
                    GMetric = GhostCheck(input_file)
                    
                    if N == "anat":
                        # Signal 2 noise ratio
                        snrCh = snrCalclualtor_chang(input_file)
                        snr_normal = snrCalclualtor_normal(input_file)   
                        
                        LMV_all = np.nan
                        GMV_all = np.nan
                        Max_mov_between_all = np.nan
                        snr_normal_vec.append(snr_normal)
                        snrCh_vec.append(snrCh)
                        
                    elif N == "diff":
                        # Signal 2 noise ratio
                        
                        snrCh = snrCalclualtor_chang(input_file)
                        snr_normal = snrCalclualtor_normal(input_file)   
                        Final,Max_mov_between,GMV,LMV = Ismotion(input_file)
                        
                        GMV_all.append(GMV)
                        LMV_all.append(LMV)
                        MI_vec_all.append(Final)
                        Max_mov_between_all.append(Max_mov_between)
                        snr_normal_vec.append(snr_normal)
                        snrCh_vec.append(snrCh)
                        
                    elif N == "func":
                        #temporal signal 2 noise ratio
                        tSNR = TsnrCalclualtor(input_file)
                        Final,Max_mov_between,GMV,LMV = Ismotion(input_file)
                        
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
            df["corresponding_img"] = img_names_new
            df['SpatRx'] = np.array(SpatRes_vec)[:,0]
            df['SpatRy'] = np.array(SpatRes_vec)[:,1]
            df['Slicethick'] = np.array(SpatRes_vec)[:,2]
            df['Ghosting'] = np.array(GMetric_vec)
            
          
            
            if N == "anat":
                 df['SNR Chang'] = np.array(snrCh_vec)
                 df['SNR Normal'] = np.array(snr_normal_vec)
                 
            elif N == "diff":
                 df['SNR Chang'] = np.array(snrCh_vec)
                 df['SNR Normal'] = np.array(snr_normal_vec)
                 df['Displacement factor (std of Mutual information)']=np.array(LMV_all)
                 #df['Maximal displacement']=AR[4]
                 
            elif N == "func":
                 df['tSNR (Averaged Brain ROI)'] = np.array(tsnr_vec)
                 df['Displacement factor (std of Mutual information)']=np.array(LMV_all)
                 #df['Maximal displacement']=AR[4]
                 
            if N=="anat":
                t2w_result= os.path.join(Path,"caculated_features_anat.csv")
                df.to_csv( t2w_result)

            elif N=="diff":    
                dti_result= os.path.join(Path,"caculated_features_diff.csv")
                df.to_csv( dti_result)   

            elif N=="func":
                fmri_result= os.path.join(Path,"caculated_features_func.csv")
                df.to_csv(fmri_result)

    if ErorrList:            
        dfNewError = pd.DataFrame(ErorrList)
        new_file_path = os.path.join(saving_path, "CanNotProcessTheseFiles.csv")
        dfNewError.to_csv(new_file_path)

    print('\n\noutput file was created:' + str(Path))
    
    print('\n\n%%%%%%%%%%%%% End of the stage 2 %%%%%%%%%%%%%%%\n\n'.upper())
    
   
#%%
