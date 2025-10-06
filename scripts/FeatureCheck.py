import os
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
import nibabel as nib
from nilearn.plotting import plot_img
from nilearn.image import index_img
import alive_progress as ap

import pv_conv2Nifti as pr
import pv_parser as par
from QC import (
    ResCalculator, GhostCheck, snrCalclualtor_chang, 
    snrCalclualtor_normal, TsnrCalclualtor, Ismotion
)


# Constants
SEQUENCE_TYPES = ['anat', 'diff', 'func']
FEATURE_COLUMNS = {
    'common': ['FileAddress', 'SpatRx', 'SpatRy', 'SpatRz', 'Ghosting'],
    'anat': ['SNR Chang', 'SNR Normal'],
    'diff': ['SNR Chang', 'SNR Normal', 'Displacement factor (std of Mutual information)'],
    'func': ['tSNR (Averaged Brain ROI)', 'Displacement factor (std of Mutual information)']
}


@dataclass
class FeatureData:
    """Container for calculated features."""
    file_paths: List[str] = field(default_factory=list)
    img_names: List[str] = field(default_factory=list)
    sequence_names: List[str] = field(default_factory=list)
    spatial_res: List[np.ndarray] = field(default_factory=list)
    ghosting: List[float] = field(default_factory=list)
    snr_chang: List[float] = field(default_factory=list)
    snr_normal: List[float] = field(default_factory=list)
    tsnr: List[float] = field(default_factory=list)
    motion_info: List[float] = field(default_factory=list)
    lmv: List[float] = field(default_factory=list)
    gmv: List[float] = field(default_factory=list)
    max_movement: List[float] = field(default_factory=list)
    
    def to_dataframe(self, seq_type: str) -> pd.DataFrame:
        """Convert feature data to pandas DataFrame."""
        df = pd.DataFrame({
            'FileAddress': self.file_paths,
            'SpatRx': [sr[0] for sr in self.spatial_res],
            'SpatRy': [sr[1] for sr in self.spatial_res],
            'SpatRz': [sr[2] for sr in self.spatial_res],
            'Ghosting': self.ghosting
        })
        
        # Add sequence-specific columns
        if seq_type == 'anat':
            if self.img_names:
                df.insert(1, 'corresponding_img', self.img_names)
            if self.sequence_names:
                df.insert(1, 'sequence name', self.sequence_names)
            df['SNR Chang'] = self.snr_chang
            df['SNR Normal'] = self.snr_normal
            
        elif seq_type == 'diff':
            if self.img_names:
                df.insert(1, 'corresponding_img', self.img_names)
            if self.sequence_names:
                df.insert(1, 'sequence name', self.sequence_names)
            df['SNR Chang'] = self.snr_chang
            df['SNR Normal'] = self.snr_normal
            df['Displacement factor (std of Mutual information)'] = self.lmv
            
        elif seq_type == 'func':
            if self.img_names:
                df.insert(1, 'corresponding_img', self.img_names)
            if self.sequence_names:
                df.insert(1, 'sequence name', self.sequence_names)
            df['tSNR (Averaged Brain ROI)'] = self.tsnr
            df['Displacement factor (std of Mutual information)'] = self.lmv
        
        return df


def load_address_files(path: str) -> Dict[str, pd.DataFrame]:
    """Load address CSV files and return dictionary by sequence type."""
    address_book = {}
    
    for file_path in Path(path).glob('*addreses*.csv'):
        filename = file_path.name
        
        if 'anat' in filename:
            address_book['anat'] = pd.read_csv(file_path)
        elif 'diff' in filename:
            address_book['diff'] = pd.read_csv(file_path)
        elif 'func' in filename:
            address_book['func'] = pd.read_csv(file_path)
    
    return address_book


def create_inspection_image(input_file: nib.Nifti1Image, 
                           save_path: str, 
                           img_name: str) -> None:
    """Create and save orthogonal view of image for manual inspection."""
    if len(input_file.shape) == 4:
        display_img = index_img(input_file, 0)
    else:
        display_img = input_file
    
    plot_img(display_img, title=img_name, output_file=save_path)


def calculate_common_features(input_file: nib.Nifti1Image) -> Tuple[np.ndarray, float]:
    """Calculate features common to all sequence types."""
    spatial_res = ResCalculator(input_file)
    ghosting = GhostCheck(input_file)
    return spatial_res, ghosting


def calculate_snr_features(input_file: nib.Nifti1Image) -> Tuple[float, float]:
    """Calculate SNR features for anatomical and diffusion sequences."""
    snr_chang = snrCalclualtor_chang(input_file)
    snr_normal = snrCalclualtor_normal(input_file)
    return snr_chang, snr_normal


def calculate_motion_features(input_file: nib.Nifti1Image) -> Tuple[float, float, float, float]:
    """Calculate motion-related features for functional and diffusion sequences."""
    return Ismotion(input_file)


def load_raw_bruker_data(file_path: str) -> Tuple[Optional[nib.Nifti1Image], Optional[str]]:
    """Load raw Bruker data and return NIfTI image and sequence key."""
    path_parts = Path(file_path).parts
    
    procno = '1'
    expno = path_parts[-1]
    study = path_parts[-2]
    raw_folder = os.sep.join(path_parts[:-2])
    proc_folder = os.path.join(raw_folder, 'proc_data')
    
    # Check for required parameter files
    visu_pars = os.path.join(file_path, 'pdata', '1', 'visu_pars')
    acqp_path = os.path.join(file_path, 'acqp')
    
    if not (os.path.isfile(visu_pars) and os.path.isfile(acqp_path)):
        raise FileNotFoundError("Missing visu_pars or acqp file")
    
    # Load Bruker data
    pv = pr.Bruker2Nifti(study, expno, procno, raw_folder, proc_folder, ftype='NIFTI_GZ')
    pv.read_2dseq(map_raw=False, pv6=False)
    input_file = nib.squeeze_image(pv.nim)
    
    # Create proper affine matrix for nilearn
    affine = np.eye(4)
    pixdim = pv.nim.header.get('pixdim')
    for i in range(4):
        affine[i, i] = pixdim[i + 1]
    
    input_file = nib.Nifti1Image(input_file.get_fdata(), affine=affine, dtype=np.int32)
    
    # Get sequence information
    scan_info = par.read_param_file(acqp_path)
    method_name = scan_info[1]["ACQ_method"].upper()
    scan_name = scan_info[1]["ACQ_scan_name"].upper()
    sequence_key = method_name + scan_name
    
    return input_file, sequence_key


def process_single_file(file_path: str, 
                       seq_type: str, 
                       output_path: str,
                       file_index: int,
                       is_raw: bool = False) -> Tuple[Optional[Dict], Optional[str]]:
    """
    Process a single file and extract features.
    
    Returns:
        Tuple of (feature_dict, error_message)
    """
    try:
        # Load data
        if is_raw:
            input_file, sequence_key = load_raw_bruker_data(file_path)
        else:
            input_file = nib.load(file_path)
            sequence_key = None
        
        # Create inspection image
        qc_path = os.path.join(output_path, "manual_slice_inspection")
        os.makedirs(qc_path, exist_ok=True)
        
        path_obj = Path(file_path)
        if is_raw:
            img_name = path_obj.parts[-2]
        else:
            img_name = f"{path_obj.parts[-2]}_{path_obj.name}"
        
        full_img_name = f"{seq_type}_{img_name}_{file_index}.png"
        full_img_name = full_img_name.replace('.nii', '').replace('.gz', '')
        
        svg_path = os.path.join(qc_path, full_img_name)
        create_inspection_image(input_file, svg_path, full_img_name)
        
        # Calculate common features
        spatial_res, ghosting = calculate_common_features(input_file)
        
        # Initialize feature dictionary
        features = {
            'file_path': file_path,
            'img_name': full_img_name,
            'sequence_key': sequence_key,
            'spatial_res': spatial_res,
            'ghosting': ghosting
        }
        
        # Calculate sequence-specific features
        if seq_type == 'anat':
            snr_chang, snr_normal = calculate_snr_features(input_file)
            features.update({
                'snr_chang': snr_chang,
                'snr_normal': snr_normal
            })
            
        elif seq_type == 'diff':
            snr_chang, snr_normal = calculate_snr_features(input_file)
            final, max_mov, gmv, lmv = calculate_motion_features(input_file)
            features.update({
                'snr_chang': snr_chang,
                'snr_normal': snr_normal,
                'motion_info': final,
                'max_movement': max_mov,
                'gmv': gmv,
                'lmv': lmv
            })
            
        elif seq_type == 'func':
            tsnr = TsnrCalclualtor(input_file)
            final, max_mov, gmv, lmv = calculate_motion_features(input_file)
            features.update({
                'tsnr': tsnr,
                'motion_info': final,
                'max_movement': max_mov,
                'gmv': gmv,
                'lmv': lmv
            })
        
        return features, None
        
    except (ValueError, SystemError, KeyError, FileNotFoundError, 
            nib.loadsave.ImageFileError) as e:
        error_msg = f"{file_path}_{type(e).__name__}"
        print(f"{type(e).__name__}: {file_path}")
        return None, error_msg


def process_sequence_type(seq_type: str, 
                          file_list: List[str], 
                          output_path: str,
                          is_raw: bool = False) -> Tuple[FeatureData, List[str]]:
    """Process all files of a given sequence type."""
    print(f'{seq_type} processing...\n')
    
    feature_data = FeatureData()
    error_list = []
    
    with ap.alive_bar(len(file_list), spinner='wait') as bar:
        for idx, file_path in enumerate(file_list, start=1):
            features, error = process_single_file(
                str(file_path), seq_type, output_path, idx, is_raw
            )
            
            if error:
                error_list.append(error)
            elif features:
                feature_data.file_paths.append(features['file_path'])
                feature_data.img_names.append(features['img_name'])
                feature_data.spatial_res.append(features['spatial_res'])
                feature_data.ghosting.append(features['ghosting'])
                
                if features.get('sequence_key'):
                    feature_data.sequence_names.append(features['sequence_key'])
                
                if 'snr_chang' in features:
                    feature_data.snr_chang.append(features['snr_chang'])
                    feature_data.snr_normal.append(features['snr_normal'])
                
                if 'tsnr' in features:
                    feature_data.tsnr.append(features['tsnr'])
                
                if 'lmv' in features:
                    feature_data.lmv.append(features['lmv'])
                    feature_data.gmv.append(features['gmv'])
                    feature_data.motion_info.append(features['motion_info'])
                    feature_data.max_movement.append(features['max_movement'])
            
            bar()
    
    return feature_data, error_list


def save_results(feature_data: FeatureData, 
                seq_type: str, 
                output_path: str) -> None:
    """Save feature data to CSV file."""
    df = feature_data.to_dataframe(seq_type)
    output_file = os.path.join(output_path, f"caculated_features_{seq_type}.csv")
    df.to_csv(output_file, index=False)


def process_features(path: str, is_raw: bool = False) -> None:
    """
    Main function to process features for all sequence types.
    
    Args:
        path: Path to the directory containing address CSV files
        is_raw: Whether processing raw Bruker data (True) or NIfTI (False)
    """
    # Load address files
    address_book = load_address_files(path)
    
    if not address_book:
        print("No address files found!")
        return
    
    all_errors = []
    
    # Process each sequence type
    for seq_type, addresses_df in address_book.items():
        if addresses_df.empty:
            continue
        
        file_list = addresses_df.iloc[:, 0].tolist()
        
        # Process files
        feature_data, errors = process_sequence_type(
            seq_type, file_list, path, is_raw
        )
        
        # Save results
        if feature_data.file_paths:
            save_results(feature_data, seq_type, path)
        
        all_errors.extend(errors)
        
        if errors:
            print(f'{len(errors)} faulty file(s) found for {seq_type}\n')
    
    # Save error list
    if all_errors:
        df_errors = pd.DataFrame({'ErrorData': all_errors})
        error_file = os.path.join(path, "CanNotProcessTheseFiles.csv")
        df_errors.to_csv(error_file, index=False)
        print(f"Error list saved to: {error_file}")
    
    print(f'\n\nOutput files created: {path}')
    print('\n\n%%%%%%%%%%%%% END OF STAGE 2 %%%%%%%%%%%%%%%\n\n')


def CheckingRawFeatures(path: str) -> None:
    """Process raw Bruker format features."""
    process_features(path, is_raw=True)


def CheckingNiftiFeatures(path: str) -> None:
    """Process NIfTI format features."""
    process_features(path, is_raw=False)
