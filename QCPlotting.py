#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 25 21:08:31 2021

@author: kalantaria
"""
import os
import pandas as pd


from openpyxl import load_workbook
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
import seaborn as sns
import matplotlib.patches as mpatches
import numpy as np
import openpyxl
import QC
#%% Path to Excel

#Path = r"\\10.209.5.114\AG_Aswendt_Projects\Student_projects\14_Aref_Kalantari_2021\Projects\CRC\QualityControl\Results\QC\BigOne_TVA_Pyramid\QuiC_Data_Result_Processed_featurs_testing.xlsx"
#Path=r"\\10.209.5.114\AG_Aswendt_Projects\Student_projects\14_Aref_Kalantari_2021\Projects\CRC\QualityControl\Datasets\Aswendt\QC\QuiC_Data_Result_Processed_featurs.xlsx"
Path = r"C:\Users\arefk\OneDrive\Desktop\BigOne_TVA_Pyramid\QuiC_Data_Result_Processed_featurs_testing.xlsx"
#%% Loading Excel table



def main():
    import argparse
    
    parser = argparse.ArgumentParser(description='Plotting QC Results')
    parser.add_argument('address',help='Address of created\
                                     Excel from CheckingFeatures_final.py (Stage II) or any other other sutibale file structure')
    args = parser.parse_args()
    Path = args.address    
    QC.QCPlot(Path)
    QC.QCtable(Path)


if __name__ == '__main__':
    main()