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
from sklearn.covariance import EllipticEnvelope
import seaborn as sns    
from sklearn.ensemble import IsolationForest
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.neighbors import LocalOutlierFactor 
from sklearn.svm import OneClassSVM
import numpy as np
from matplotlib.pyplot import imshow
import os
import pandas as pd
import glob
from openpyxl import load_workbook
import matplotlib.backends.backend_pdf
import seaborn as sns
import matplotlib.patches as mpatches
import numpy as np
import openpyxl
import nibabel as nib
import glob
import openpyxl
import nibabel as nii
import os
import time
from nibabel.testing import data_path
import matplotlib.pyplot as plt
#from skimage import data
#from skimage.filters import try_all_threshold
#from skimage.filters import threshold_isodata
from scipy import ndimage
from scipy import signal
import matplotlib.pyplot as plt
from matplotlib.transforms import Affine2D
from matplotlib import transforms
import mpl_toolkits.axisartist.floating_axes as floating_axes
import changSNR as ch
import pv_conv2Nifti as pr
import alive_progress as ap
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
#%% Goasing 
def GoastCheck(input_file):
    
    #input_file= nib.load(tf)
    img = input_file
    img_data = img.get_fdata()
    img_shape = np.shape(img_data)
    MI_vec = []
    n = 1
    Mmos = []
    while (img_shape[1]%(2**n)) == 0:
        Mmos.append(img_shape[1]/2**n)
        n = n+1
    Mmos = np.asarray(Mmos)
    if len(img_shape)>3:
        img_data = np.mean(img_data,axis=-1)
                
    
    Im_ref = img_data[:,:,int(img_shape[2]/2)]
    for ii in range(0,int(img_shape[1])):
        Im_rol = np.roll(Im_ref,ii)
        MI_vec.append(mutualInfo(Im_rol,Im_ref))
        
    peaks_strong, prop = signal.find_peaks(MI_vec, height = 0.25*max(MI_vec))
    peaks_weak, prop = signal.find_peaks(MI_vec)
    
    StrongGoast = np.sum(np.isin(peaks_strong,Mmos))
    WeekGoast = np.sum(np.isin(peaks_weak,Mmos))
    
    if WeekGoast > 2 or StrongGoast > 0:
        GMetric = "yes"
    else:
        GMetric = "no"
    
    
    #plt.plot((MI_vec))
    #plt.show()
    
    return GMetric
  

#%% Res function

def Image_Selection(input_file):
    
    img = input_file
    img_data=img.get_fdata()
    img_shape=img_data.shape 
    middle=int(img_shape[2]/2)
    if len(img_shape) == 3:
            selected_img= img_data[:, :,middle]
            
    elif len(img_shape) > 3:
        time_middle=int(img_shape[3]/2)
        selected_img= img_data[:, :,middle,time_middle]
        
    selected_img= np.rot90(selected_img,k=-1)
    return selected_img


#%% 


def ResCalculator(input_file):
    
    HDR = input_file.header
    Spati_Res = HDR['pixdim'][1:4]
    
    return Spati_Res



#%% SNR function
def snrCalclualtor_chang(input_file):

    imgData = input_file
    IM = np.asanyarray(imgData.dataobj)
    imgData = np.squeeze(np.ndarray.astype(IM, 'float64'))
    
    
    snr_chang_slice_vec = []
    ns = imgData.shape[2]  # Number of slices
    n_dir = imgData.shape[-1]  # Number of directions if dti dataset
    nd = imgData.ndim
    
    ns_lower = int(np.floor(ns/2) - 2)
    ns_upper = int(np.floor(ns/2) + 2)

    #print('/NewData/',end=" ")
    #for slc in range(ns_lower,ns_upper):
    for slc in range(ns_lower,ns_upper):    
        #   Print % of progress
        #print('S' + str(slc + 1), end=",")

        # Decision if the input data is DTI type or T2w
        if nd == 3:
            Slice = imgData[:, :, slc]
            try:
                curSnrCHMap, estStdChang, estStdChangNorm = ch.calcSNR(Slice, 0, 1)
            except ValueError:
                estStdChang = np.nan
            snr_chang_slice = 20 * np.log10(np.mean(Slice)/estStdChang)
            snr_chang_slice_vec.append(snr_chang_slice)
        else:
            for bb in range(5,n_dir):
                Slice = imgData[:, :,slc,bb]
                try:
                    curSnrCHMap, estStdChang, estStdChangNorm = ch.calcSNR(Slice, 0, 1)
                except ValueError:
                    estStdChang = np.nan
                snr_chang_slice = 20 * np.log10(np.mean(Slice)/estStdChang)
                snr_chang_slice_vec.append(snr_chang_slice)
        
    snr_chang_slice_vec = np.array(snr_chang_slice_vec)    
    snrCh = np.mean(snr_chang_slice_vec[~np.isinf(snr_chang_slice_vec) *  ~np.isnan(snr_chang_slice_vec)])

    return snrCh

#%% SNR function 2
def snrCalclualtor_normal(input_file):
    
    IM = np.asanyarray(input_file.dataobj)
    imgData = np.squeeze(np.ndarray.astype(IM, 'float64'))
    Data = imgData
    
    S = np.shape(np.squeeze(Data))
    if len(S) == 3:
        imgData = np.squeeze(Data)
    if len(S) == 4:
        imgData = np.squeeze(Data[:,:,:,int((S[2]/2))])
    
    S = np.shape(np.squeeze(imgData))
    
    #local thresholding
    imgData_new = np.zeros(S[0:3]);
# =============================================================================
#     for ii in range(0,S[2]):
#         temp_image = imgData[:,:,ii]
#         global_thresh = threshold_isodata(temp_image)
#         binary_global = temp_image > global_thresh
#         imgData_new[:,:,ii] = binary_global
        
# =============================================================================
    
    COM=[int(i) for i in (ndimage.measurements.center_of_mass(imgData))]
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
   for i, Slice in enumerate(slices):
       axes[i].imshow(Slice.T, cmap="gray", origin="lower")
       

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

#%% TSNR function

def TsnrCalclualtor(input_file):
    imgData = input_file
    IM = np.asanyarray(imgData.dataobj)
    imgData = np.ndarray.astype(IM, 'float64')
    signal_averge_over_time = imgData[:,:,:,10:].mean(axis=-1) 
    signal_std_over_time = imgData[:,:,:,10:].std(axis=-1) 
    tSNR_map = 20 * np.log10(signal_averge_over_time/signal_std_over_time)
    
    S = np.shape(input_file)
     #local thresholding
    imgData_new = np.zeros(S[0:3])
    imgData_average = np.mean(imgData,axis=-1)
# =============================================================================
#     for ii in range(0,S[2]):
#         temp_image = imgData_average[:,:,ii]
#         global_thresh = threshold_isodata(temp_image)
#         binary_global = temp_image > global_thresh
#         imgData_new[:,:,ii] = binary_global
#         
# =============================================================================
    
    COM=[int(i) for i in (ndimage.measurements.center_of_mass(imgData_average))]
    r = np.floor(0.10*(np.mean([S[0:2]])))
    Mask = sphere(S[0:3], int(r) , COM)
    tSNR = np.mean(tSNR_map[Mask])
    
    return tSNR


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
    GMV=[]
    imgData = input_file
    IM = np.asanyarray(imgData.dataobj)
    imgData = np.ndarray.astype(IM[:,:,:,10:], 'float64')
    S = np.shape(imgData)
    temp_mean = imgData.mean(axis=(0,1,3))
    temp_max = temp_mean.argmax()
    temp_Data = imgData[:,:,temp_max,:]
    Im_fix = temp_Data[:,:,0]
    Im_rot = temp_Data
    
    MI_all = []
    for z in range(1,S[-1]):
        
        MI = mutualInfo(Im_fix,Im_rot[:,:,z])
        MI_all.append(MI)
    
    Final = np.asarray(MI_all)
    Max_mov_between = str([Final.argmin()+10,Final.argmax()+10])
    GMV = getrange(Final)
    LMV = np.std(Final)
    
    return Final,Max_mov_between,GMV,LMV

#%% Getting range 

def getrange(numbers):
    return max(numbers) - min(numbers)


#%% Plotting QC Histogram and etc.

def QCPlot(Path):
    
    saving_path = (Path) 
    QC_fig_path =os.path.join( (Path) , "QCfigures")
    if not os.path.isdir(QC_fig_path):
        os.mkdir(QC_fig_path)

       
    Abook = []
    Names =[]
    for file in glob.glob(os.path.join(Path, '*caculated_features*.csv')) :
        
        if "DTI" in file.upper():
            dti_path= file
            dti_features= pd.read_csv(dti_path)
            Abook.append(dti_features)
            Names.append("DTI")
        elif "FMRI" in file.upper():
            fmri_path= file
            fmri_features= pd.read_csv(fmri_path)
            Abook.append(fmri_features)
            Names.append("rsfMRI")
        elif "T2W" in file.upper():    
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
    hh = 1
    rr = 1
    for nn,N in enumerate(Names):

            
        d= Abook[nn]
        COL = Abook[nn].columns
        
        for cc,C in enumerate(COL):
           
            Data = list(d[C])
            
            
            if C == 'SNR Chang' or C == 'tSNR (Averaged Brain ROI)' or C =='SNR Normal':
             
             #plot histogrm
                plt.figure(hh,figsize=(10,5))
                ax2 = plt.subplot(1,1,1)
                for dd,DD in enumerate(Data):
                    if DD == np.inf:
                        Data[dd] = np.nan
                
                q75, q25 = np.nanpercentile(Data, [75 ,25])
                iqr = q75 - q25
                
                B = round((np.nanmax(Data)-np.nanmin(Data)) / (2 * iqr / (len(Data)**(1/3))))
                if B*5 > 25:
                    BB = 25
                else:
                    BB = 25
                
                y, x, bars = plt.hist(Data, bins= B*7, histtype= 'bar',edgecolor='white')
                plt.xlabel(N+': '+C + ' [a.u.]')
                plt.ylabel("Frequency")
                ax2.spines['right'].set_visible(False)
                ax2.spines['top'].set_visible(False)
                plt.locator_params(axis='x', nbins=BB)
                #calculate interquartile range of values in the 'points' column
                
                if C == 'Displacement factor (std of Mutual information)':
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
                
                
    
                red_patch = mpatches.Patch(color='red', label='Discard')
                blue_patch = mpatches.Patch(color='tab:blue', label='Keep')
                plt.legend(handles=[blue_patch,red_patch])
                plt.savefig(os.path.join(QC_fig_path,C+N+".tiff"),dpi=250)
                plt.close()
                   # plt.savefig(os.path.dirname(Path) + "\ResHomogenity.png",dpi=300)
               
    hh = hh + 1           
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
    plt.savefig(os.path.join(QC_fig_path,"Spatial_Resolution.tiff"),dpi=250)
    plt.close()
    

#%% Adjusting the existing feature table by adding a new sheet to it with the data that need to be discarded

def QCtable(Path):
    #Path= r"C:\Users\Erfan\Downloads\Documents"
    ML_algorythms= ML(Path)
    ML_algorythms=pd.concat(ML_algorythms) 
    ML_algorythms[['One_class_SVM',' EllipticEnvelope','IsolationForest',"LocalOutlierFactor"]]=ML_algorythms[['One_class_SVM',' EllipticEnvelope','IsolationForest',"LocalOutlierFactor"]]==-1 
    Abook = []
    Names =[]
    for file in glob.glob(os.path.join(Path, '*caculated_features*.csv')) :
        
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
            
            
            if C == 'SNR Chang' or C == 'tSNR (Averaged Brain ROI)' or C =='SNR Normal':
                
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
                
                
            if C == 'Displacement factor (std of Mutual information)':
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
 

         
         

    
 
    
    #prepare ML outliers

    ml_outliers= ML_algorythms[(ML_algorythms["One_class_SVM" ]==True) | (ML_algorythms["IsolationForest" ]==True) |
                     (ML_algorythms["LocalOutlierFactor" ]==True) | (ML_algorythms[" EllipticEnvelope" ]==True)  ]
 
    
    
    ml_outliers= ml_outliers.rename(columns={"address": "Pathes"})  
    ml_outliers["Problematic Quality Feature"  ]= len(ml_outliers)*["ml_outlier" ] 
    ml_outliers["Sequence Type"  ]= len(ml_outliers)*["unknown" ]     
    ml_outliers= ml_outliers[["Pathes","Sequence Type","Problematic Quality Feature",'One_class_SVM',' EllipticEnvelope',
                              'IsolationForest',"LocalOutlierFactor"]]
                            
    
    Overlap=[]
    for ml_outlier in ML_algorythms["address"]:
        for classic_outlier in Pathes :
            
            if ml_outlier== classic_outlier:
                Overlap.append(ml_outlier)
    for path in Overlap :
        Overlap=ML_algorythms[ML_algorythms['address'] == path]
        
        
    Overlap= Overlap.rename(columns={"address": "Pathes"})  
     
    List = {"Pathes":Pathes,"Sequence Type":ST, "Problematic Quality Feature":COE}
    
    df = pd.DataFrame(List)
    
    
    df =df.merge(Overlap[['Pathes','One_class_SVM',' EllipticEnvelope','IsolationForest',"LocalOutlierFactor"]])
    df =df.append(ml_outliers)
    

    ML_number=list(df[["One_class_SVM" ,'IsolationForest',"LocalOutlierFactor",' EllipticEnvelope']].sum(axis=1))
    
    ML_number= [True if x>=3 else False for x in ML_number]
            
        
    
    df["final_outliers"]=ML_number
    final_result = os.path.join(Path,"unaccountable_data.csv")
    df.to_csv( final_result)    
 
    
 
    
 
#%%  
def ML(Path) :

    #prepare data
    #Path= r"C:\Users\Erfan\Downloads\Documents"
 
    
    csv_path=["caculated_features_DTI.csv","caculated_features_T2w.csv","caculated_features_fMRI.csv"]
    result=[]
    for N, csv in enumerate(glob.glob(os.path.join(Path, '*_features_*.csv'))):
        
       csv_path=os.path.join(Path,csv)
       Abook= pd.read_csv(csv_path)
       Abook= Abook.dropna()
       address= [i for i in Abook.iloc[:,1]]
       X = [[i] for i in  Abook.iloc[:,6]]
       
############## Fit the One-Class SVM 
       nu = 0.05
       gamma = 2.0
       clf = OneClassSVM(gamma="auto", kernel="poly", nu=nu,shrinking=False).fit(X)
       svm_pre =clf.predict(X)
############## EllipticEnvelope
        
       elpenv = EllipticEnvelope(contamination=0.025, random_state=1)
       ell_pred = elpenv.fit_predict(X)
    
############## IsolationForest
   
       iforest = IsolationForest(n_estimators=100, max_samples='auto', 
                              contamination=0.05, max_features=1.0, 
                              bootstrap=False, n_jobs=-1, random_state=1)
       iso_pred = iforest.fit_predict(X)
    
############## LocalOutlierFactor
    
       lof = LocalOutlierFactor(n_neighbors=20, algorithm='auto',
                             metric='minkowski', contamination=0.04,
                             novelty=False, n_jobs=-1)
       local_pred = lof.fit_predict(X)
    
        
############## saving result
       algorythms=[svm_pre,ell_pred,iso_pred,local_pred]
       result.append(algorythms)
       result[N]= np.dstack((result[N][0], result[N][1],result[N][2],result[N][3]))
       result[N]= result[N][0]
       result[N]= pd.DataFrame(result[N], columns = ['One_class_SVM',' EllipticEnvelope','IsolationForest',"LocalOutlierFactor"])
       
       result[N]["address"] = address
        
    return(result)


#%% For Questions please Contact: aref.kalantari-sarcheshmeh@uk-koeln.de