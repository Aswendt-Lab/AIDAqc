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
from pathlib import Path
from typing import Dict, List, Tuple
import argparse

import pandas as pd
import alive_progress as ap
import pv_parser as par
import QC
import FeatureCheck as fc


# Constants
DTI_KEYWORDS = ["DTI", "STRUCT", "DWI", "DIFFUS"]
FMRI_KEYWORDS = ["RESTING", "FUN", "RSF", "FMRI", "BOLD", "RS-"]
T2_KEYWORDS = ["T2W", "T1W", "ANAT", "RARE", "TURBO", "T1_F", "T2_F"]
NOT_ALLOWED = ["LOC", "PIL", "FISP", "WOB", "NOIS", "SINGL", "MRS", "B0M", "FIELD"]
TYPE_STRINGS = ['diff', 'func', 'anat']


def print_header():
    """Print welcome message and contact information."""
    print("Hello! Are you ready to get rid of bad quality data?")
    print('------------------------------------------------------------')
    print('Thank you for using our code. Contact:')
    print('aref.kalantari-sarcheshmeh@uk-koeln.de / markus.aswendt@uk-koeln.de')
    print('Lab: AG Neuroimaging and Neuroengineering of Experimental Stroke, University Hospital Cologne')
    print('Web: https://neurologie.uk-koeln.de/forschung/ag-neuroimaging-neuroengineering/')
    print('------------------------------------------------------------')


def print_footer():
    """Print closing message."""
    print('------------------------------------------------------------')
    print('Thank you for using our code. For questions please contact us via:')
    print('aref.kalantari-sarcheshmeh@uk-koeln.de or markus.aswendt@uk-koeln.de')
    print('Lab: AG Neuroimaging and neuroengineering of experimental stroke University Hospital Cologne')
    print('Web:https://neurologie.uk-koeln.de/forschung/ag-neuroimaging-neuroengineering/')
    print('------------------------------------------------------------')


def initialize_address_book() -> Dict[str, List[str]]:
    """Initialize the address book dictionary."""
    return {type_str: [] for type_str in TYPE_STRINGS}


def check_keywords(text: str, keywords: List[str]) -> bool:
    """Check if any keyword is present in the text."""
    return any(keyword in text for keyword in keywords)


def classify_sequence_type(key: str, param_data: Dict) -> Tuple[str, bool]:
    """
    Classify the sequence type based on keywords.
    
    Returns:
        Tuple of (sequence_type, is_allowed)
    """
    flag_anat = check_keywords(key, T2_KEYWORDS)
    flag_struct = check_keywords(key, DTI_KEYWORDS)
    flag_func = check_keywords(key, FMRI_KEYWORDS)
    flag_not_allowed = check_keywords(key, NOT_ALLOWED)
    flag_epi = "EPI" in key
    
    if flag_not_allowed:
        return None, False
    
    if flag_struct:
        return "diff", True
    elif flag_func:
        return "func", True
    elif flag_anat and not flag_epi:
        return "anat", True
    elif flag_epi:
        # Additional check for EPI sequences
        time_points = param_data.get("ACQ_time_points", [])
        movie_frames = param_data.get("ACQ_n_movie_frames", 0)
        if movie_frames != len(time_points):
            return "func", True
        else:
            return "diff", True
    
    return None, False


def parse_raw_files(initial_path: str) -> Tuple[List[str], int]:
    """Parse raw MR files and return list of acqp files."""
    path_pattern = os.path.join(initial_path, "**", "acqp")
    
    with ap.alive_bar(title='Parsing through folders...', length=10, 
                      stats=False, monitor=False) as bar:
        text_files = glob.glob(path_pattern, recursive=True)
        bar()
    
    print(f'TOTAL NUMBER OF {len(text_files)} FILES WERE FOUND: PARSING FINISHED!')
    return text_files, len(text_files)


def extract_raw_sequences(text_files: List[str]) -> Tuple[Dict[str, List[str]], List[str]]:
    """Extract usable sequences from raw files."""
    address_book = initialize_address_book()
    error_list = []
    seen_dates = set()
    count = 0
    
    with ap.alive_bar(len(text_files), 
                      title='EXTRACTING T1/T2 WEIGHTED, DIFF AND FUNC SEQUENCES:',
                      length=10, stats=False, spinner='wait') as bar:
        
        for filepath in text_files:
            try:
                scan_info, param_data = par.read_param_file(filepath)
                method_name = param_data["ACQ_method"].upper()
                protocol_name = param_data["ACQ_protocol_name"].upper()
                key = method_name + protocol_name
                scan_date = scan_info['Date']
                
            except (KeyError, UnicodeDecodeError) as e:
                print(f"{type(e).__name__}: {filepath}")
                error_list.append(filepath)
                bar()
                continue
            
            # Skip duplicates based on date
            if scan_date in seen_dates:
                bar()
                continue
            
            seq_type, is_allowed = classify_sequence_type(key, param_data)
            
            if seq_type and is_allowed:
                address_book[seq_type].append(os.path.dirname(filepath))
                count += 1
            
            seen_dates.add(scan_date)
            bar()
    
    duplicates_removed = len(text_files) - len(seen_dates)
    print(f' {count} FILES WERE EXTRACTED!')
    print(f' {duplicates_removed} DUPLICATES WERE ELIMINATED!')
    
    return address_book, error_list


def parse_nifti_files(initial_path: str, suffix: str) -> List[str]:
    """Parse NIfTI files and return list of file paths."""
    patterns = [
        os.path.join(initial_path, "**", f"*{suffix}.nii.gz"),
        os.path.join(initial_path, "**", f"*{suffix}.nii")
    ]
    
    with ap.alive_bar(title='Parsing through folders...', length=10,
                      stats=False, monitor=False) as bar:
        text_files = []
        for pattern in patterns:
            text_files.extend(glob.glob(pattern, recursive=True))
        bar()
    
    print(f'TOTAL NUMBER OF {len(text_files)} FILES WERE FOUND: PARSING FINISHED!')
    return text_files


def extract_nifti_sequences(text_files: List[str]) -> Dict[str, List[str]]:
    """Extract usable sequences from NIfTI files."""
    address_book = initialize_address_book()
    
    for filepath in text_files:
        path_parts = filepath.split(os.sep)
        
        # Check path components from filename backwards
        seq_type = None
        for idx in range(len(path_parts)):
            part_upper = path_parts[-idx - 1].upper()
            
            flag_anat = check_keywords(part_upper, T2_KEYWORDS)
            flag_struct = check_keywords(part_upper, DTI_KEYWORDS)
            flag_func = check_keywords(part_upper, FMRI_KEYWORDS)
            flag_not_allowed = check_keywords(part_upper, NOT_ALLOWED)
            
            if any([flag_anat, flag_struct, flag_func, flag_not_allowed]):
                break
        
        # Determine sequence type
        if flag_not_allowed:
            continue
        
        if flag_struct:
            seq_type = "diff"
        elif flag_func:
            seq_type = "func"
        elif flag_anat:
            seq_type = "anat"
        
        if seq_type:
            address_book[seq_type].append(filepath)
    
    return address_book


def save_address_book(address_book: Dict[str, List[str]], 
                      saving_path: str, 
                      format_type: str):
    """Save address book to CSV files."""
    for type_str, addresses in address_book.items():
        if addresses:
            df = pd.DataFrame(addresses, columns=[0])
            filename = f"{format_type}_data_addreses_{type_str}.csv"
            filepath = os.path.join(saving_path, filename)
            df.to_csv(filepath, sep=',', index=False)


def save_error_list(error_list: List[str], saving_path: str):
    """Save error list to CSV file."""
    if error_list:
        df_error = pd.DataFrame({'ErrorData': error_list})
        error_path = os.path.join(saving_path, "CanNotOpenTheseFiles.csv")
        df_error.to_csv(error_path, index=False)
        print("Some data could not be opened - Check CanNotOpenTheseFiles.csv")


def cleanup_and_organize(saving_path: str):
    """Clean up temporary files and organize results."""
    # Remove address files
    for file in glob.glob(os.path.join(saving_path, '*data_addreses*.csv')):
        os.remove(file)
    
    # Create and move calculated features
    features_dir = os.path.join(saving_path, "calculated_features")
    os.makedirs(features_dir, exist_ok=True)
    
    for old_file in glob.glob(os.path.join(saving_path, '*caculated_features*.csv')):
        filename = os.path.basename(old_file)
        new_file = os.path.join(features_dir, filename)
        os.replace(old_file, new_file)


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Parser of all MR files: This code will parse through every "
                   "possible folder behind a defined initial path, looking for MR "
                   "data files of any type. Then it will extract the wanted files "
                   "and eliminate any duplicates."
    )
    
    parser.add_argument('-i', '--initial_path', required=True,
                       help='Initial path to start the parsing')
    parser.add_argument('-o', '--output_path', required=True,
                       help='Path where the results should be saved')
    parser.add_argument('-f', '--format_type', required=True,
                       choices=["nifti", "raw"],
                       help='Format of your dataset: nifti or raw Bruker')
    parser.add_argument('-s', '--suffix', default="",
                       help='Suffix for data files (e.g., test for test.nii.gz)')
    parser.add_argument('-e', '--exclude', nargs='+',
                       help='Sequences to exclude from analysis')
    
    return parser.parse_args()


def main():
    """Main execution function."""
    args = parse_arguments()
    
    # Update NOT_ALLOWED list with exclusions
    if args.exclude:
        NOT_ALLOWED.extend([e.upper() for e in args.exclude])
    
    # Setup
    QC.tic()
    print_header()
    os.makedirs(args.output_path, exist_ok=True)
    
    # Parse files based on format
    if args.format_type == "raw":
        text_files, _ = parse_raw_files(args.initial_path)
        address_book, error_list = extract_raw_sequences(text_files)
        save_error_list(error_list, args.output_path)
    else:  # nifti
        text_files = parse_nifti_files(args.initial_path, args.suffix)
        address_book = extract_nifti_sequences(text_files)
    
    # Save results
    save_address_book(address_book, args.output_path, args.format_type)
    print(f'\n\nCSV files were created: {args.output_path}')
    print('\n\n%%%%%%%%%%%%% END OF STAGE 1 %%%%%%%%%%%%%%%')
    
    # Stage 2: Feature calculation
    print('\nSTARTING STAGE 2...')
    print('\nCALCULATING FEATURES...\n')
    print('This will take some time depending on the dataset size.\n\n')
    
    if args.format_type == "raw":
        fc.CheckingRawFeatures(args.output_path)
    else:
        fc.CheckingNiftiFeatures(args.output_path)
    
    QC.toc()
    
    # Plotting and cleanup
    print('PLOTTING QUALITY FEATURES...\n')
    QC.QCPlot(args.output_path)
    QC.QCtable(args.output_path, args.format_type)
    
    cleanup_and_organize(args.output_path)
    
    print('\n\n%%%%%%%%%%%%%QUALITY FEATURE PLOTS SUCCESSFULLY CREATED%%%%%%%%%%%%%%%\n\n')
    print_footer()


if __name__ == "__main__":
    main()
