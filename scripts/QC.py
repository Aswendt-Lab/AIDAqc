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
from sklearn.ensemble import IsolationForest
from sklearn.neighbors import LocalOutlierFactor 
from sklearn.svm import OneClassSVM
import numpy as np
import os
import pandas as pd
import glob
from datetime import datetime
import subprocess
import hashlib
import matplotlib.patches as mpatches
import time
import matplotlib.pyplot as plt
from scipy import ndimage
from scipy import signal
import changSNR as ch
from matplotlib.ticker import MaxNLocator
from matplotlib import font_manager as fm
from matplotlib.font_manager import FontProperties
import re
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
#%% Ghosting 
def GhostCheck(input_file):
    
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
    
    StrongGhost = np.sum(np.isin(peaks_strong,Mmos))
    WeekGhost = np.sum(np.isin(peaks_weak,Mmos))
    
    if WeekGhost > 2 or StrongGhost > 0:
        GMetric = True
    else:
        GMetric = False
    
    
    #plt.plot((MI_vec))
    #plt.show()
    
    return GMetric


#%% Ghosting numeric metric (for MI_diff.csv / MI_func.csv extra outputs)
def GhostCheckValue(input_file):
    """Return a numeric ghosting score derived from the same peak-count logic as GhostCheck.

    - GhostCheck() returns a boolean (ghosting present / not present).
    - This function returns an integer score:
        score = WeekGhost + StrongGhost

    Where WeekGhost/StrongGhost are the counts of MI peaks that fall on the expected
    ghost offsets (Mmos). This keeps the *decision logic* unchanged while exposing
    a magnitude-like value for reporting.
    """

    img = input_file
    img_data = img.get_fdata()
    img_shape = np.shape(img_data)
    MI_vec = []

    n = 1
    Mmos = []
    while (img_shape[1] % (2 ** n)) == 0:
        Mmos.append(img_shape[1] / 2 ** n)
        n = n + 1
    Mmos = np.asarray(Mmos)

    if len(img_shape) > 3:
        img_data = np.mean(img_data, axis=-1)

    Im_ref = img_data[:, :, int(img_shape[2] / 2)]
    for ii in range(0, int(img_shape[1])):
        Im_rol = np.roll(Im_ref, ii)
        MI_vec.append(mutualInfo(Im_rol, Im_ref))

    peaks_strong, _ = signal.find_peaks(MI_vec, height=0.25 * max(MI_vec))
    peaks_weak, _ = signal.find_peaks(MI_vec)

    StrongGhost = int(np.sum(np.isin(peaks_strong, Mmos)))
    WeekGhost = int(np.sum(np.isin(peaks_weak, Mmos)))

    return WeekGhost + StrongGhost
  



#%% Ghosting intensity-based metric (Ghost-to-Signal Ratio, GSR)
def GhostGSRValue(input_file, percentile=70, ghost_shift_fraction=0.5, phase_axis=1):
    """Return a continuous ghosting score using an intensity-based Ghost-to-Signal Ratio (GSR).

    This is designed for reporting in the MI_* extra outputs only (does not affect the main pipeline logic).

    Method (single mid-slice):
      1) Create a foreground/object mask using a percentile threshold on positive intensities.
      2) Shift that mask by ~half FOV along the phase-encode axis to where EPI-like ghosts appear.
      3) Measure mean background-corrected absolute intensity in the shifted (ghost-only) region
         relative to the mean object intensity.

    Returns:
      GhostingScore = mean(|I_ghost - median(I_bg)|) / (mean(I_obj) + eps)

    Notes:
      - For 4D data, uses the temporal mean image (same as GhostCheck).
      - phase_axis defaults to 1 to match the original GhostCheck loop over img_shape[1].
    """

    img = input_file
    img_data = img.get_fdata()
    img_shape = np.shape(img_data)

    if len(img_shape) > 3:
        img_data = np.mean(img_data, axis=-1)

    # Mid-slice
    slc = img_data[:, :, int(img_shape[2] / 2)]
    slc = np.asarray(slc, dtype=np.float64)

    pos = slc[slc > 0]
    if pos.size < 10:
        return np.nan

    thr = np.percentile(pos, percentile)
    obj_mask = slc > thr

    # Fallback if the mask is too small
    if np.sum(obj_mask) < 10:
        thr = np.percentile(pos, 50)
        obj_mask = slc > thr

    if np.sum(obj_mask) < 10:
        return np.nan

    # Determine shift (half FOV by default)
    try:
        n_phase = slc.shape[phase_axis]
    except Exception:
        return np.nan

    shift = int(round(n_phase * ghost_shift_fraction))
    if shift == 0:
        shift = n_phase // 2

    ghost_mask = np.roll(obj_mask, shift=shift, axis=phase_axis)
    ghost_only = ghost_mask & (~obj_mask)

    if np.sum(ghost_only) < 10:
        # No valid ghost-only region
        return 0.0

    bg = slc[~obj_mask]
    bg_med = float(np.median(bg)) if bg.size else 0.0

    s = float(np.mean(slc[obj_mask]))
    g = float(np.mean(np.abs(slc[ghost_only] - bg_med)))

    eps = 1e-12
    return g / (s + eps)
#%% Res function


def ResCalculator(input_file):
    
    HDR = input_file.header
    Spati_Res = HDR['pixdim'][1:4]
    
    return Spati_Res



#%% SNR function
def snrCalclualtor_chang(input_file):

    imgData = input_file
    IM = np.asanyarray(imgData.dataobj)
    imgData = np.squeeze(np.ndarray.astype(IM, 'float64'))

    mm = imgData.mean()
    if mm == 0:
        snrCh = np.nan
        return snrCh
    """
    cc = 0
    while mm < 1:
        mm = mm *10
        cc = cc+1
    imgData = imgData * (10**cc)
    """
    
    Sone = len(imgData.shape)
    if Sone < 3:
        snrCh = np.nan
        return snrCh
        
    
    snr_chang_slice_vec = []
    ns = imgData.shape[2]  # Number of slices
    n_dir = imgData.shape[-1]  # Number of directions if dti dataset
    if len(imgData.shape) > 3:    
        if n_dir < 10 :
            fff = 0
            print()
            print("Warning: Be aware that the size of the 4th dimension (difusion direction or timepoints) is less than 10. This might result in unstable values")
            print()
        else:
            fff = 5
        
    nd = imgData.ndim
    if ns > 4:
        ns_lower = int(np.floor(ns/2) - 2)
        ns_upper = int(np.floor(ns/2) + 2)
    else:
        ns_lower=0
        ns_upper=ns
        
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
            for bb in range(fff,n_dir-1):
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
def snrCalclualtor_normal(
    input_file,
    output_dir=None,
    save_sphere_png=True,
    sphere_png_name=None,
    sphere_radius_scale=1.0,
    sphere_radius_factor=0.5,
    use_ellipsoid_if_needed=False,
    ellipsoid_z_scale=1.0,
    ellipsoid_xy_factor=0.25,
):
    
    
    
    IM = np.asanyarray(input_file.dataobj)
    imgData = np.squeeze(np.ndarray.astype(IM, 'float64'))
    
    Sone = len(imgData.shape)
    if Sone < 3:
        imgData = np.tile(imgData[:, :, np.newaxis], (1, 1, 10))
    
    Data = imgData
    
    S = np.shape(np.squeeze(Data))
    #print(S)
    if len(S) == 3:
        imgData = np.squeeze(Data)
    if len(S) == 4:
        imgData = np.squeeze(Data[:,:,:,0]) #int((S[-1]/2))
    
    S = np.shape(np.squeeze(imgData))
    
    #local thresholding
    #imgData_new = np.zeros(S[0:3]);
# =============================================================================
#     for ii in range(0,S[2]):
#         temp_image = imgData[:,:,ii]
#         global_thresh = threshold_isodata(temp_image)
#         binary_global = temp_image > global_thresh
#         imgData_new[:,:,ii] = binary_global
        
# =============================================================================
    
    COM=[int(i) for i in (ndimage.measurements.center_of_mass(imgData))]
    # Sphere radius in voxels (dynamic, based on in-plane size).
    # For thin stacks (small Z), a true 3D sphere is limited by Z.
    # If use_ellipsoid_if_needed=True, we keep a large in-plane radius and
    # limit only the Z semisize (ellipsoid/spheroid).
    base_dim = float(np.mean(S[0:2]))
    r_xy = np.floor(float(sphere_radius_factor) * base_dim * float(sphere_radius_scale))
    r_xy = int(max(1, r_xy))

    # maximum possible *sphere* radius that fits in all dimensions
    max_r_sphere = max(1, int(np.floor(min(S) / 2)))
    r_sphere = int(min(r_xy, max_r_sphere))

    if use_ellipsoid_if_needed and r_sphere < r_xy:
        # Thin anatomical volume: use an ellipsoid so that the in-plane ROI
        # is not limited by the small Z dimension. Keep this ROI smaller than
        # the previous 0.5 * in-plane dimension fallback.
        r_xy_ellipsoid = int(max(
            1,
            np.floor(float(ellipsoid_xy_factor) * base_dim * float(sphere_radius_scale))
        ))
        max_rz = max(1, int(np.floor(S[2] / 2)))
        rz = int(np.ceil(r_xy_ellipsoid * float(ellipsoid_z_scale)))
        rz = int(max(1, min(max_rz, rz)))
        Mask = sphere(
            S,
            r_xy_ellipsoid,
            COM,
            semisizes=(r_xy_ellipsoid, r_xy_ellipsoid, rz)
        )
        r_used_for_meta = r_xy_ellipsoid
    else:
        Mask = sphere(S, r_sphere, COM)
        r_used_for_meta = r_sphere
    Singal = np.mean(imgData[Mask])
    
    
    x = int(np.ceil(S[0]*0.15))
    y = int(np.ceil(S[1]*0.15))
    z = int(np.ceil(S[2]*0.15))
    
    MaskN = np.zeros(S[0:3]);
    MaskN[:x,:y,:z] = 2
    MaskN[:x,-y:,:z] = 2
    MaskN[-x:,:y,:z] = 2
    MaskN[-x:,-y:,:z] = 2
    MaskN[:x,:y,-z:] = 2
    MaskN[:x,-y:,-z:] = 2
    MaskN[-x:,:y,-z:] = 2
    MaskN[-x:,-y:,-z:] = 2


    # Save the sphere mask + overlay (including edge/corner noise ROIs) to manual_slice_inspection
    # Enabled by default. If output_dir is not provided, it is inferred from the input file path.
    if save_sphere_png:
        output_dir = _resolve_output_dir(output_dir, input_file)

        # Use the same naming logic as manual slice inspection: derive from the input filename
        if sphere_png_name is None:
            stem = _infer_image_stem(input_file)
            sphere_png_name = f"{stem}_sphere_mask_snr_normal.png"

        base, ext = os.path.splitext(sphere_png_name)
        overlay_name = (base.replace('mask', 'overlay') + ext) if ext else (base.replace('mask', 'overlay') + '.png')

        try:
            sphere_dir = output_dir
            os.makedirs(sphere_dir, exist_ok=True)

            # avoid overwriting when multiple scans are processed
            sphere_png_name = _make_unique_filename(sphere_dir, sphere_png_name)
            base, ext = os.path.splitext(sphere_png_name)
            overlay_name = (base.replace('mask', 'overlay') + ext) if ext else (base.replace('mask', 'overlay') + '.png')
            overlay_name = _make_unique_filename(sphere_dir, overlay_name)

            # Edge/corner ROIs used for noise estimation
            edge_bool = (MaskN > 0)

            save_sphere_mask_png(Mask, sphere_dir, filename=sphere_png_name, center=COM, radius=int(r_used_for_meta))
            # Additional overlay on anatomical middle slices (orthogonal views)
            save_sphere_overlay_png(
                imgData,
                Mask,
                sphere_dir,
                filename=overlay_name,
                edge_mask=edge_bool
            )
        except Exception as e:
            print(f'Warning: could not save sphere mask/overlay PNG: {e}')
    
    
    
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
    SNR = 20 * np.log10(Singal/Noise_std)
    if np.isinf(SNR):
        SNR = np.nan
        print("Impossible: Infinite values were the result of SNR")
        print("Possible reason: already ROI extracted/preprocessed data with zeros around the ROI. S/0=inf'")
        print("for continuity, inf is replaced with NaN ...")
    return SNR



def show_slices(slices):
   """ Function to display row of image slices """
   fig, axes = plt.subplots(1, len(slices))
   for i, Slice in enumerate(slices):
       axes[i].imshow(Slice.T, cmap="gray", origin="lower")
       

def sphere(shape, radius, position, semisizes=None):
    """Generate an n-dimensional spherical mask."""
    # assume shape and position have the same length and contain ints
    # the units are pixels / voxels (px for short)
    # radius is a int or float in px
    assert len(position) == len(shape)
    n = len(shape)
    if semisizes is None:
        semisizes = (radius,) * len(shape)
    else:
        assert len(semisizes) == len(shape)

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




def _infer_image_stem(input_file):
    """Infer a stable stem/name for the current image for file naming.

    Tries nibabel-style get_filename(); falls back to 'image'.
    Returns a filesystem-safe stem (no extension).
    """
    if QC_DEFAULT_NAME_PREFIX:
        return re.sub(r"[^A-Za-z0-9._-]+", "_", str(QC_DEFAULT_NAME_PREFIX)).strip("._-") or "image"

    path = None
    try:
        if hasattr(input_file, 'get_filename') and callable(getattr(input_file, 'get_filename')):
            path = input_file.get_filename()
    except Exception:
        path = None

    if not path:
        # Some nibabel objects may store filename differently
        for attr in ('filename', '_filename'):
            try:
                v = getattr(input_file, attr, None)
                if isinstance(v, str) and v:
                    path = v
                    break
            except Exception:
                pass

    if not path:
        stem = 'image'
    else:
        base = os.path.basename(str(path))
        # handle .nii.gz explicitly
        if base.lower().endswith('.nii.gz'):
            stem = base[:-7]
        else:
            stem = os.path.splitext(base)[0]

    # make filesystem safe
    stem = re.sub(r'[^A-Za-z0-9._-]+', '_', stem).strip('._-')
    return stem or 'image'


def _infer_output_dir(input_file):
    """Infer an output directory from the input file.

    If the image is loaded from disk (nibabel), we use the directory containing the NIfTI.
    Returns None if it cannot be inferred.
    """
    path = None
    try:
        if hasattr(input_file, 'get_filename') and callable(getattr(input_file, 'get_filename')):
            path = input_file.get_filename()
    except Exception:
        path = None

    if not path:
        # Some nibabel objects may store filename differently
        for attr in ('filename', '_filename'):
            try:
                v = getattr(input_file, attr, None)
                if isinstance(v, str) and v:
                    path = v
                    break
            except Exception:
                pass

    if not path:
        return None
    try:
        d = os.path.dirname(str(path))
        return d if d else None
    except Exception:
        return None


def _make_unique_filename(folder, filename):
    """If filename exists in folder, append _02, _03, ..."""
    folder = str(folder)
    base, ext = os.path.splitext(filename)
    candidate = filename
    i = 2
    while os.path.exists(os.path.join(folder, candidate)):
        candidate = f"{base}_{i:02d}{ext}"
        i += 1
    return candidate


# --- Default output directory handling (for sphere PNG exports) ---
# If your pipeline writes calculated_features*.csv into a specific output folder,
# call set_qc_output_dir(<that folder>) ONCE before processing to make sure all
# sphere PNG exports go there as well.
QC_DEFAULT_OUTPUT_DIR = None
QC_DEFAULT_NAME_PREFIX = None  # optional, set by pipeline for consistent per-scan naming

def set_qc_name_prefix(prefix):
    """Set a filename prefix used for sphere PNG outputs (without extension).

    This is useful when input_file has no on-disk filename (e.g., created in-memory),
    and you want the sphere PNGs to match the manual_slice_inspection naming scheme.
    """
    global QC_DEFAULT_NAME_PREFIX
    QC_DEFAULT_NAME_PREFIX = str(prefix) if prefix else None


def set_qc_output_dir(path):
    """Set a default output directory used for QC PNG exports."""
    global QC_DEFAULT_OUTPUT_DIR
    QC_DEFAULT_OUTPUT_DIR = str(path) if path is not None else None


def _resolve_output_dir(output_dir, input_file):
    """Resolve output directory for PNG exports.

    Priority:
      1) explicit output_dir argument
      2) QC_DEFAULT_OUTPUT_DIR set via set_qc_output_dir(...)
      3) infer from input_file filename (directory containing NIfTI)
      4) current working directory (last resort)
    """
    if output_dir is not None:
        return str(output_dir)
    if QC_DEFAULT_OUTPUT_DIR is not None:
        return str(QC_DEFAULT_OUTPUT_DIR)
    inferred = _infer_output_dir(input_file)
    if inferred is not None:
        return str(inferred)
    return os.getcwd()


def _manual_slice_dir(base_dir):
    """Return manual_slice_inspection directory under base_dir (or base_dir if already that folder)."""
    base_dir = str(base_dir)
    if os.path.basename(os.path.normpath(base_dir)) == 'manual_slice_inspection':
        return base_dir
    return os.path.join(base_dir, 'manual_slice_inspection')

def save_sphere_mask_png(mask, output_dir, filename='sphere_mask.png', center=None, radius=None, dpi=200):
    """Save a 3D spherical boolean mask as a PNG (3 orthogonal slices).

    Parameters
    ----------
    mask : ndarray (3D, bool or 0/1)
        The spherical mask returned by `sphere(...)`.
    output_dir : str or path-like
        Folder where the PNG will be saved (e.g., the same folder as calculated_features*.csv).
    filename : str
        Output PNG filename.
    center : tuple of 3 ints, optional
        Sphere center (x, y, z). If None, uses the middle of the volume.
    radius : int or float, optional
        Sphere radius (voxels). Only used for the figure title.
    dpi : int
        PNG resolution.

    Returns
    -------
    str
        Full path to the written PNG file.
    """
    if output_dir is None:
        raise ValueError('output_dir must be provided to save the sphere mask PNG.')

    mask = np.asarray(mask)
    if mask.ndim != 3:
        raise ValueError(f'mask must be 3D, got shape {mask.shape}')

    os.makedirs(output_dir, exist_ok=True)

    # Choose a center slice for display
    if center is None:
        center = tuple(int(s // 2) for s in mask.shape)
    cx, cy, cz = (int(center[0]), int(center[1]), int(center[2]))
    cx = max(0, min(cx, mask.shape[0] - 1))
    cy = max(0, min(cy, mask.shape[1] - 1))
    cz = max(0, min(cz, mask.shape[2] - 1))

    # Convert to uint8 for clean rendering
    msk = mask.astype(np.uint8)

    fig, axes = plt.subplots(1, 3, figsize=(9, 3), dpi=dpi)
    axes[0].imshow(msk[:, :, cz].T, cmap='gray', origin='lower')
    axes[0].set_title(f'Axial (z={cz})', fontsize=9)

    axes[1].imshow(msk[:, cy, :].T, cmap='gray', origin='lower')
    axes[1].set_title(f'Coronal (y={cy})', fontsize=9)

    axes[2].imshow(msk[cx, :, :].T, cmap='gray', origin='lower')
    axes[2].set_title(f'Sagittal (x={cx})', fontsize=9)

    for ax in axes:
        ax.axis('off')

    title = 'Sphere mask'
    if radius is not None:
        title += f' (r={radius} vox)'
    fig.suptitle(title, fontsize=10)
    fig.tight_layout()

    out_path = os.path.join(output_dir, filename)
    fig.savefig(out_path, bbox_inches='tight')
    plt.close(fig)
    return out_path


def save_sphere_overlay_png(anat_image, mask, output_dir, filename='sphere_overlay.png', dpi=200, edge_mask=None):
    """Save an overlay PNG (orthogonal anatomical views + ROI contours).

    Requested behavior: show the *middle slices* of the anatomical volume.
    Practical behavior: if the ROI does not intersect a middle slice in a given view,
    we fall back to the closest slice that *does* intersect the ROI (so the contour
    is visible and you can verify placement).

    Overlays:
      - `mask` (central ROI) as a solid contour on the selected slice.
      - `edge_mask` (edge/corner noise ROIs) as a dashed contour.

    Note on `edge_mask` visualization:
      The noise ROIs live in the *corners* of the volume and often only exist
      in the first/last few slices. If we draw them only on the same slice
      index used for the central ROI, they frequently won't intersect that
      slice and the contour will be empty (especially for thin-Z anat stacks).
      To make them reliably visible, we plot edge ROIs using a projection
      (any-voxel) along the viewing axis.
    """

    if output_dir is None:
        raise ValueError('output_dir must be provided to save the sphere overlay PNG.')

    anat = np.asarray(anat_image)
    msk = np.asarray(mask).astype(bool)
    if anat.ndim != 3:
        raise ValueError(f'anat_image must be 3D, got shape {anat.shape}')
    if msk.ndim != 3:
        raise ValueError(f'mask must be 3D, got shape {msk.shape}')
    if anat.shape != msk.shape:
        raise ValueError(f'anat_image and mask shapes must match. Got {anat.shape} vs {msk.shape}')

    edge = None
    if edge_mask is not None:
        edge = np.asarray(edge_mask).astype(bool)
        if edge.ndim != 3:
            raise ValueError(f'edge_mask must be 3D, got shape {edge.shape}')
        if edge.shape != anat.shape:
            raise ValueError(f'edge_mask shape must match anat_image. Got {edge.shape} vs {anat.shape}')

    os.makedirs(output_dir, exist_ok=True)

    # Preferred (middle) slice indices
    mid_x, mid_y, mid_z = (anat.shape[0] // 2, anat.shape[1] // 2, anat.shape[2] // 2)

    def _choose_slice(mask3d, axis, preferred):
        """Prefer `preferred`, else choose closest slice with any mask voxels."""
        if axis == 0:
            proj = mask3d.any(axis=(1, 2))
        elif axis == 1:
            proj = mask3d.any(axis=(0, 2))
        else:
            proj = mask3d.any(axis=(0, 1))

        preferred = int(max(0, min(preferred, len(proj) - 1)))
        if proj[preferred]:
            return preferred, True
        idxs = np.where(proj)[0]
        if idxs.size == 0:
            return preferred, False
        best = int(idxs[np.argmin(np.abs(idxs - preferred))])
        return best, False

    # Choose slice indices (prefer middle; fall back if needed)
    mx, used_mid_x = _choose_slice(msk, axis=0, preferred=mid_x)
    my, used_mid_y = _choose_slice(msk, axis=1, preferred=mid_y)
    mz, used_mid_z = _choose_slice(msk, axis=2, preferred=mid_z)

    # Prepare slices (transpose for display like your other viewers)
    ax_img = anat[:, :, mz].T
    ax_msk = msk[:, :, mz].T
    co_img = anat[:, my, :].T
    co_msk = msk[:, my, :].T
    sa_img = anat[mx, :, :].T
    sa_msk = msk[mx, :, :].T

    if edge is not None:
        # Use projections so corner/cubicle ROIs are visible even if they
        # don't intersect the selected slice index.
        # Axial view shows X-Y plane -> project over Z
        ax_edge = edge.any(axis=2).T
        # Coronal view shows X-Z plane -> project over Y
        co_edge = edge.any(axis=1).T
        # Sagittal view shows Y-Z plane -> project over X
        sa_edge = edge.any(axis=0).T
    else:
        ax_edge = co_edge = sa_edge = None

    fig, axes = plt.subplots(1, 3, figsize=(9, 3), dpi=dpi)

    # Axial
    axes[0].imshow(ax_img, cmap='gray', origin='lower')
    if ax_msk.any():
        axes[0].contour(ax_msk.astype(float), levels=[0.5], colors='r', linewidths=1.0)
    if ax_edge is not None and ax_edge.any():
        axes[0].contour(ax_edge.astype(float), levels=[0.5], colors='y', linewidths=1.2, linestyles='--')
    axes[0].set_title(f'Axial z={mz}' + (' (mid)' if used_mid_z else ' (nearest ROI)'), fontsize=9)

    # Coronal
    axes[1].imshow(co_img, cmap='gray', origin='lower')
    if co_msk.any():
        axes[1].contour(co_msk.astype(float), levels=[0.5], colors='r', linewidths=1.0)
    if co_edge is not None and co_edge.any():
        axes[1].contour(co_edge.astype(float), levels=[0.5], colors='y', linewidths=1.2, linestyles='--')
    axes[1].set_title(f'Coronal y={my}' + (' (mid)' if used_mid_y else ' (nearest ROI)'), fontsize=9)

    # Sagittal
    axes[2].imshow(sa_img, cmap='gray', origin='lower')
    if sa_msk.any():
        axes[2].contour(sa_msk.astype(float), levels=[0.5], colors='r', linewidths=1.0)
    if sa_edge is not None and sa_edge.any():
        axes[2].contour(sa_edge.astype(float), levels=[0.5], colors='y', linewidths=1.2, linestyles='--')
    axes[2].set_title(f'Sagittal x={mx}' + (' (mid)' if used_mid_x else ' (nearest ROI)'), fontsize=9)

    for ax in axes:
        ax.axis('off')

    fig.suptitle('Sphere overlay (mid slices with ROI fallback)', fontsize=10)
    fig.tight_layout()

    out_path = os.path.join(output_dir, filename)
    fig.savefig(out_path, bbox_inches='tight')
    plt.close(fig)
    return out_path

#%% TSNR function

def TsnrCalclualtor(
    input_file,
    output_dir=None,
    save_sphere_png=True,
    sphere_png_name=None,
    sphere_radius_scale=1.0,
    sphere_radius_factor=0.50,
    use_ellipsoid_if_needed=False,
    ellipsoid_z_scale=1.0,
):
    imgData = input_file
    IM = np.asanyarray(imgData.dataobj)
    S=IM.shape
    if len(S) == 3:    
        IM = IM.reshape((S[0],S[1],1,S[2]))
    
    
    imgData = np.ndarray.astype(IM, 'float64')
    if IM.shape[-1] < 10:
        fff = 0
    else:
        fff = 10
 
    signal_averge_over_time = imgData[:,:,:,fff:].mean(axis=-1) 
    signal_std_over_time = imgData[:,:,:,fff:].std(axis=-1) 
    tSNR_map = 20 * np.log10(signal_averge_over_time/signal_std_over_time)
    
    S = np.shape(IM)
     #local thresholding
    #imgData_new = np.zeros(S[0:3])
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
    # Sphere radius in voxels (dynamic, based on in-plane size).
    base_dim = float(np.mean(S[0:2]))
    r_xy = np.floor(float(sphere_radius_factor) * base_dim * float(sphere_radius_scale))
    r_xy = int(max(1, r_xy))

    max_r_sphere = max(1, int(np.floor(min(S[0:3]) / 2)))
    r_sphere = int(min(r_xy, max_r_sphere))

    if use_ellipsoid_if_needed and r_sphere < r_xy:
        max_rz = max(1, int(np.floor(S[2] / 2)))
        rz = int(np.ceil(r_xy * float(ellipsoid_z_scale)))
        rz = int(max(1, min(max_rz, rz)))
        Mask = sphere(S[0:3], r_xy, COM, semisizes=(r_xy, r_xy, rz))
        r_used_for_meta = r_xy
    else:
        Mask = sphere(S[0:3], r_sphere, COM)
        r_used_for_meta = r_sphere

    # Save the sphere mask + overlay to manual_slice_inspection
    # Enabled by default. If output_dir is not provided, it is inferred from the input file path.
    if save_sphere_png:
        output_dir = _resolve_output_dir(output_dir, input_file)

        # Use the same naming logic as manual slice inspection: derive from the input filename
        if sphere_png_name is None:
            stem = _infer_image_stem(input_file)
            sphere_png_name = f"{stem}_sphere_mask_tsnr.png"
        base, ext = os.path.splitext(sphere_png_name)
        overlay_name = (base.replace('mask', 'overlay') + ext) if ext else (base.replace('mask', 'overlay') + '.png')
        try:
            sphere_dir = output_dir
            os.makedirs(sphere_dir, exist_ok=True)

            # avoid overwriting when multiple scans are processed
            sphere_png_name = _make_unique_filename(sphere_dir, sphere_png_name)
            base, ext = os.path.splitext(sphere_png_name)
            overlay_name = (base.replace('mask', 'overlay') + ext) if ext else (base.replace('mask', 'overlay') + '.png')
            overlay_name = _make_unique_filename(sphere_dir, overlay_name)

            save_sphere_mask_png(Mask, sphere_dir, filename=sphere_png_name, center=COM, radius=int(r_used_for_meta))
            # Additional overlay on anatomical middle slices (orthogonal views)
            save_sphere_overlay_png(
                imgData_average,
                Mask,
                sphere_dir,
                filename=overlay_name
            )
        except Exception as e:
            print(f'Warning: could not save sphere mask PNG: {e}')
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


#%% Motion detection of rsFRI function (based on mutual information)

def Ismotion(input_file):
    GMV=[]
    imgData = input_file
    IM = np.asanyarray(imgData.dataobj)
    S = IM.shape
    if len(S) == 3:    
        IM = IM.reshape((S[0],S[1],1,S[2]))
    
    
    if IM.shape[-1] < 11 :
        fff = 0
    else:
        fff = 10
 
    imgData = np.ndarray.astype(IM[:,:,:,fff:], 'float64')
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
    QC_fig_path = os.path.join( (Path) , "QCfigures")
    if not os.path.isdir(QC_fig_path):
        os.mkdir(QC_fig_path)

       
    Abook = []
    Names =[]
    for file in glob.glob(os.path.join(Path, '*calculated_features*.csv')) :
        
        if "diff" in file:
            dti_path= file
            dti_features= pd.read_csv(dti_path)
            Abook.append(dti_features)
            Names.append("diff")
        elif "func" in file:
            fmri_path= file
            fmri_features= pd.read_csv(fmri_path)
            Abook.append(fmri_features)
            Names.append("func")
        elif "anat" in file:    
             t2w_path= file
             t2w_features= pd.read_csv(t2w_path)
             Abook.append(t2w_features)
             Names.append("anat")    

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
    # Set font properties
    title_font = {'family': 'serif', 'fontname': 'DejaVu Sans'}
    label_font = {'family': 'serif', 'fontname': 'DejaVu Sans'}
    tick_font = {'family': 'serif', 'fontname': 'DejaVu Sans'}
    
    for nn, N in enumerate(Names):
        COL = list(Abook[nn].columns)
        COL.pop(0)
        D = Abook[nn]
        
        for cc, C in enumerate(COL):
            Data = list(D[C])
            
            if C == 'SNR Chang' or C == 'tSNR (Averaged Brain ROI)' or C == 'SNR Normal' or C == 'Displacement factor (std of Mutual information)':
                # Plot histogram
                cm = 1/2.54  # centimeters in inches
                plt.figure(hh, figsize=(9, 5), dpi=300)
                ax2 = plt.subplot(1, 1, 1, label="histogram")
                
                for dd, DD in enumerate(Data):  # If tSNR and SNR chang are also adjusted, this section can be eliminated
                    if DD == np.inf:
                        Data[dd] = np.nan
                
                q75, q25 = np.nanpercentile(Data, [75, 25])
                iqr = q75 - q25
                
                B = round((np.nanmax(Data) - np.nanmin(Data)) / (2 * iqr / (len(Data) ** (1/3))))
                if B * 5 > 22:
                    XX = 22
                else:
                    XX = B * 5
                
                y, x, bars = plt.hist(Data, bins=B * 7, histtype='bar', edgecolor='white')
                plt.xlabel(N + ': ' + C + ' [a.u.]', fontdict=label_font)
                plt.ylabel("Frequency", fontdict=label_font)
                ax2.spines['right'].set_visible(False)
                ax2.spines['top'].set_visible(False)
                plt.locator_params(axis='x', nbins=XX)
                ax2.yaxis.set_major_locator(MaxNLocator(integer=True))
                
                # Calculate interquartile range of values in the 'points' column
                if C == 'Displacement factor (std of Mutual information)':
                    ll = q75 + 1.5 * iqr
                    plt.text(1.07 * ll, 2 * max(y) / 3, 'Q3 + 1.5*IQ', color='grey', fontdict=label_font)
                    for b, bar in enumerate(bars):
                        if bar.get_x() > ll:
                            bar.set_facecolor("red")
                else:
                    ll = q25 - 1.5 * iqr
                    plt.text(1.001 * ll, 2 * max(y) / 3, 'Q1 - 1.5*IQ', color='grey', fontdict=label_font)
                    for b, bar in enumerate(bars):
                        if bar.get_x() < ll:
                            bar.set_facecolor("red")
                
                plt.axvline(ll, color='grey', linestyle='--')
                plt.suptitle(N + ': ' + C, fontdict=title_font)
                
                red_patch = mpatches.Patch(color='red', label='Discard')
                blue_patch = mpatches.Patch(color='tab:blue', label='Keep')
                # Modify the legend with smaller font size and Times New Roman font
                legend = plt.legend(handles=[blue_patch, red_patch], fontsize=8)
                #legend.get_frame().set_linewidth(0.0)  # Remove legend border

                # Set Times New Roman font for legend text
                # Set Times New Roman font for legend text
                for text in legend.get_texts():
                   text.set_fontfamily('serif')
                   text.set_fontsize(8)

                
                # Set the font for axis ticks
                ax2.xaxis.set_tick_params(labelsize=8)
                ax2.yaxis.set_tick_params(labelsize=8)
                
                base_filename = os.path.join(QC_fig_path, C + N)
                plt.savefig(base_filename + ".png", dpi=300)
                plt.savefig(base_filename + ".svg", format='svg')

                plt.close()
                
        hh = hh + 1
    
    plt.figure(hh, figsize=(9, 5), dpi=300)
    for nn, N in enumerate(Names):
        COL = list(Abook[nn].columns)
        COL.pop(0)
        D = Abook[nn]
        for cc, C in enumerate(COL):
            Data = list(D[C])
            if C == 'SpatRx' or C == 'SpatRy' or C == 'SpatRz':
                # Plot pie plots
                labels = list(set(Data))
                sizes = [Data.count(l) for l in labels]
                labels = list(np.round(labels, 3))
                labels2 = [str(l) + ' mm' for l in labels]
                
                ax1 = plt.subplot(len(Names), 3, rr)
                ax1.pie(sizes, labels=labels2, autopct='%1.0f%%', startangle=180)
                ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
                ax1.set_title(N + ':' + C, fontdict=title_font)
                plt.suptitle('Resolution homogeneity between data', weight="bold")
                
                # Set the font for axis ticks
                ax1.xaxis.set_tick_params(labelsize=8)
                ax1.yaxis.set_tick_params(labelsize=8)
                
                rr = rr + 1
    
    base_filename = os.path.join(QC_fig_path, "Spatial_Resolution")
    plt.savefig(base_filename + ".png", dpi=300)
    plt.savefig(base_filename + ".svg", format='svg')
    plt.close()

#%%
# machine learning methods
def ML(Path, format_type) :

    result=[]
    for N, csv in enumerate(glob.glob(os.path.join(Path, '*_features_*.csv'))):
        csv_path = csv
        csv_path=os.path.join(Path,csv)
        Abook= pd.read_csv(csv_path)
        if np.any(Abook.isnull().all()[:]):
            print("The following csv file contains NaN values for one or more of its features:")
            print(csv_path)
            print("Voting can not be conducted.")
            print("Analyzing next sequence...")
            continue
                 
        Abook= Abook.dropna(how='all',axis='columns')
        Abook= Abook.dropna(how='any')
        address= [i for i in Abook.iloc[:,1]]
        if format_type == "raw":
            sequence_name = [i for i in Abook.iloc[:,2]]
            img_name = [i for i in Abook.iloc[:,3]]
            X =  Abook.iloc[:,7:]
        elif format_type == "nifti":
            img_name = [i for i in Abook.iloc[:,2]]
            X =  Abook.iloc[:,6:]

       #X=preprocessing.normalize(X)
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
        if "diff" in csv:
            dti=["diff"]*len(result[N])
            result[N]["sequence_type"] = dti
           
        elif "func" in csv:
            fmri=["func"]*len(result[N])
            result[N]["sequence_type"] = fmri          
       
        elif "anat" in csv :
            t2w=["anat"]*len(result[N])
            result[N]["sequence_type"] = t2w    
       
        result[N]["Pathes"] = address
        if format_type == "raw":
            result[N]["sequence_name"] = sequence_name
        result[N]["corresponding_img"] = img_name
        
    return(result)


#%% Adjusting the existing feature table by adding a new sheet to it with the data that need to be discarded

def _qa_sha256(file_path):
    """Return SHA256 for provenance without changing the source file."""
    h = hashlib.sha256()
    with open(file_path, "rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def _qa_git_commit(start_path):
    """Return the current Git commit when the output resides in/under a repository."""
    candidates = [os.path.abspath(start_path), os.path.dirname(os.path.abspath(__file__))]
    for candidate in candidates:
        try:
            return subprocess.check_output(
                ["git", "-C", candidate, "rev-parse", "--short", "HEAD"],
                stderr=subprocess.DEVNULL,
                text=True,
            ).strip()
        except Exception:
            pass
    return "not available"


def _qa_modality_from_name(file_path):
    name = os.path.basename(file_path).lower()
    if "anat" in name:
        return "anat"
    if "diff" in name:
        return "diff"
    if "func" in name:
        return "func"
    return "unknown"


def _qa_find_id_column(df):
    """Prefer the explicit AIDAqc address column, otherwise use the second CSV column."""
    preferred = ["FileAddress", "Pathes", "File Address", "fileaddress"]
    for col in preferred:
        if col in df.columns:
            return col
    if len(df.columns) >= 2:
        return df.columns[1]
    return df.columns[0]


def _qa_numeric(value, digits=4):
    try:
        if pd.isna(value):
            return "-"
        return f"{float(value):.{digits}f}"
    except Exception:
        return str(value)


def GenerateQAReport(Path, format_type):
    """
    Create a self-contained AIDAqc MRI QA protocol after Stage III.

    The report reads existing calculated_features_*.csv, MI_*.csv, votings.csv,
    and QCfigures/*.png files. It does not recalculate features or change any
    QC/outlier decisions.
    """
    from reportlab.lib import colors
    from reportlab.lib.enums import TA_CENTER, TA_LEFT
    from reportlab.lib.pagesizes import A4
    from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
    from reportlab.lib.units import mm
    from reportlab.platypus import (
        SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, PageBreak,
        Image as RLImage, KeepTogether
    )

    output_path = os.path.abspath(Path)
    report_path = os.path.join(output_path, "AIDAqc_QA_Report.pdf")

    # ---------- Read the already-generated Stage I/II/III outputs ----------
    feature_files = sorted(glob.glob(os.path.join(output_path, "calculated_features_*.csv")))
    mi_files = sorted(glob.glob(os.path.join(output_path, "MI_*.csv")))
    voting_path = os.path.join(output_path, "votings.csv")
    voting_df = pd.read_csv(voting_path) if os.path.isfile(voting_path) else pd.DataFrame()

    feature_tables = {}
    for file_path in feature_files:
        modality = _qa_modality_from_name(file_path)
        if modality != "unknown":
            feature_tables[modality] = pd.read_csv(file_path)

    mi_tables = {}
    for file_path in mi_files:
        modality = _qa_modality_from_name(file_path)
        if modality != "unknown":
            mi_tables[modality] = pd.read_csv(file_path)

    # ---------- Styles ----------
    styles = getSampleStyleSheet()
    title_style = ParagraphStyle(
        "AIDAqcTitle", parent=styles["Title"], fontName="Helvetica-Bold",
        fontSize=20, leading=24, alignment=TA_CENTER, spaceAfter=8
    )
    subtitle_style = ParagraphStyle(
        "AIDAqcSubtitle", parent=styles["Normal"], fontName="Helvetica",
        fontSize=10, leading=13, alignment=TA_CENTER, textColor=colors.HexColor("#444444")
    )
    h1 = ParagraphStyle(
        "AIDAqcH1", parent=styles["Heading1"], fontName="Helvetica-Bold",
        fontSize=14, leading=17, spaceBefore=5, spaceAfter=8
    )
    h2 = ParagraphStyle(
        "AIDAqcH2", parent=styles["Heading2"], fontName="Helvetica-Bold",
        fontSize=11, leading=14, spaceBefore=4, spaceAfter=5
    )
    body = ParagraphStyle(
        "AIDAqcBody", parent=styles["BodyText"], fontName="Helvetica",
        fontSize=8.5, leading=11.5, spaceAfter=5
    )
    small = ParagraphStyle(
        "AIDAqcSmall", parent=body, fontSize=7, leading=9, textColor=colors.HexColor("#555555")
    )
    scan_title = ParagraphStyle(
        "AIDAqcScanTitle", parent=h1, fontSize=13, leading=16
    )

    def _footer(canvas, doc):
        canvas.saveState()
        canvas.setFont("Helvetica", 7)
        canvas.setFillColor(colors.HexColor("#666666"))
        canvas.drawString(18 * mm, 10 * mm, "AIDAqc MRI Quality Assurance Protocol")
        canvas.drawRightString(192 * mm, 10 * mm, f"Page {doc.page}")
        canvas.restoreState()

    doc = SimpleDocTemplate(
        report_path, pagesize=A4,
        rightMargin=16 * mm, leftMargin=16 * mm,
        topMargin=16 * mm, bottomMargin=16 * mm,
        title="AIDAqc MRI Quality Assurance Protocol",
        author="AIDAqc"
    )
    story = []

    # ---------- Cover / run summary ----------
    story.append(Spacer(1, 18 * mm))
    story.append(Paragraph("AIDAqc", title_style))
    story.append(Paragraph("MRI Quality Assurance Protocol", title_style))
    story.append(Paragraph(
        "Automated report generated after Stage III outlier detection. "
        "All values shown here originate from outputs already produced by AIDAqc; "
        "report generation does not modify the QC analysis.", subtitle_style
    ))
    story.append(Spacer(1, 10 * mm))

    generated = datetime.now().astimezone().strftime("%Y-%m-%d %H:%M:%S %Z")
    git_commit = _qa_git_commit(output_path)
    run_rows = [
        ["Output directory", output_path],
        ["Input format", str(format_type)],
        ["Report generated", generated],
        ["AIDAqc Git commit", git_commit],
        ["Feature tables", str(len(feature_files))],
        ["Stage III result", "votings.csv present" if os.path.isfile(voting_path) else "votings.csv missing"],
    ]
    run_table = Table(run_rows, colWidths=[42 * mm, 125 * mm], repeatRows=0)
    run_table.setStyle(TableStyle([
        ("FONTNAME", (0, 0), (-1, -1), "Helvetica"),
        ("FONTSIZE", (0, 0), (-1, -1), 8),
        ("BACKGROUND", (0, 0), (0, -1), colors.HexColor("#EEEEEE")),
        ("FONTNAME", (0, 0), (0, -1), "Helvetica-Bold"),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("GRID", (0, 0), (-1, -1), 0.35, colors.HexColor("#CCCCCC")),
        ("LEFTPADDING", (0, 0), (-1, -1), 5),
        ("RIGHTPADDING", (0, 0), (-1, -1), 5),
        ("TOPPADDING", (0, 0), (-1, -1), 4),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
    ]))
    story.append(run_table)
    story.append(Spacer(1, 8 * mm))

    summary_rows = [["Modality", "Scans", "Stage III flagged", "Flagged (%)"]]
    total_scans = 0
    total_flagged = 0
    for modality in ["anat", "diff", "func"]:
        df = feature_tables.get(modality)
        n = 0 if df is None else len(df)
        total_scans += n
        flagged = 0
        if not voting_df.empty and "sequence_type" in voting_df.columns:
            flagged = int((voting_df["sequence_type"].astype(str).str.lower() == modality).sum())
        total_flagged += flagged
        pct = (100.0 * flagged / n) if n else 0.0
        summary_rows.append([modality, n, flagged, f"{pct:.1f}"])
    summary_rows.append(["Total", total_scans, total_flagged,
                         f"{(100.0 * total_flagged / total_scans):.1f}" if total_scans else "0.0"])
    summary_table = Table(summary_rows, colWidths=[35 * mm, 35 * mm, 48 * mm, 40 * mm], repeatRows=1)
    summary_table.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#D9EAF7")),
        ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
        ("FONTNAME", (0, -1), (-1, -1), "Helvetica-Bold"),
        ("FONTSIZE", (0, 0), (-1, -1), 8),
        ("ALIGN", (1, 1), (-1, -1), "CENTER"),
        ("GRID", (0, 0), (-1, -1), 0.35, colors.HexColor("#BDBDBD")),
        ("TOPPADDING", (0, 0), (-1, -1), 5),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
    ]))
    story.append(Paragraph("Stage III summary", h2))
    story.append(summary_table)
    story.append(Paragraph(
        "A scan is counted as Stage III flagged when it appears in votings.csv, i.e. at least one "
        "of the five Stage III outlier indicators was positive. No additional pass/fail threshold is "
        "introduced by this report.", small
    ))
    story.append(PageBreak())

    # ---------- Methods / protocol ----------
    story.append(Paragraph("QA methods and decision logic", h1))
    story.append(Paragraph(
        "AIDAqc performs feature-based MRI quality assessment followed by Stage III outlier detection. "
        "The report summarizes the calculated feature tables, the modality-specific MI/GSR outputs when "
        "available, and the final Stage III voting table.", body
    ))

    methods_rows = [
        ["Measure", "Use in AIDAqc / report"],
        ["SNR", "Signal-to-noise features from the calculated feature tables. Low-value statistical outliers are identified using the existing IQR rule."],
        ["tSNR", "Temporal SNR where available. Low-value statistical outliers are identified using the existing IQR rule."],
        ["Motion", "For 4D diffusion/functional data, motion variability is represented by the standard deviation of the mutual-information time series."],
        ["Ghosting", "The main calculated feature table retains the existing AIDAqc ghosting result."],
        ["GSR", "The additional MI_*.csv files provide a continuous intensity-based ghost-to-signal ratio for descriptive QA analysis."],
        ["Stage III", "One-class SVM, Elliptic Envelope, Isolation Forest, Local Outlier Factor, and the feature-based statistical outlier indicator are combined as five independent votes."],
    ]
    methods_rows_wrapped = [[Paragraph(str(c), small) for c in r] for r in methods_rows]
    methods_table = Table(methods_rows_wrapped, colWidths=[34 * mm, 132 * mm], repeatRows=1)
    methods_table.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#D9EAF7")),
        ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
        ("FONTNAME", (0, 1), (0, -1), "Helvetica-Bold"),
        ("FONTSIZE", (0, 0), (-1, -1), 7.5),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("GRID", (0, 0), (-1, -1), 0.35, colors.HexColor("#C8C8C8")),
        ("TOPPADDING", (0, 0), (-1, -1), 5),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
    ]))
    story.append(methods_table)
    story.append(Spacer(1, 5 * mm))
    story.append(Paragraph("Stage III statistical outlier rules", h2))
    story.append(Paragraph(
        "For SNR Chang, SNR Normal, and tSNR, the existing Stage III statistical method flags values below "
        "Q1 - 1.5 x IQR. For the displacement factor (standard deviation of mutual information), values above "
        "Q3 + 1.5 x IQR are flagged. These rules are reported as implemented; the PDF generator does not "
        "change or recompute them.", body
    ))
    story.append(PageBreak())

    # ---------- Dataset-level QC figures ----------
    story.append(Paragraph("Dataset-level QA figures", h1))
    figure_files = sorted(glob.glob(os.path.join(output_path, "QCfigures", "*.png")))
    if figure_files:
        for fig_path in figure_files:
            try:
                img = RLImage(fig_path)
                max_w, max_h = 170 * mm, 100 * mm
                scale = min(max_w / img.imageWidth, max_h / img.imageHeight, 1.0)
                img.drawWidth = img.imageWidth * scale
                img.drawHeight = img.imageHeight * scale
                caption = os.path.splitext(os.path.basename(fig_path))[0].replace("_", " ")
                story.append(KeepTogether([
                    Paragraph(caption, h2),
                    img,
                    Spacer(1, 4 * mm)
                ]))
            except Exception as fig_error:
                story.append(Paragraph(f"Figure could not be embedded: {os.path.basename(fig_path)} ({fig_error})", small))
    else:
        story.append(Paragraph(
            "No PNG files were found in QCfigures at report-generation time. The quantitative report remains complete; "
            "run QCPlot before Stage III if aggregate QC figures should be embedded.", body
        ))
    story.append(PageBreak())

    # ---------- ROI placement overlays ----------
    story.append(Paragraph("ROI placement overlays", h1))
    story.append(Paragraph(
        "The following images show the automatically positioned signal ROI "
        "and, where applicable, the peripheral background/noise ROIs used for "
        "the SNR calculation. These figures are generated during feature "
        "calculation and are included here for visual verification of ROI placement.",
        body
    ))

    roi_overlay_files = sorted(set(
        glob.glob(os.path.join(output_path, "*_sphere_overlay_*.png")) +
        glob.glob(os.path.join(output_path, "manual_slice_inspection", "*_sphere_overlay_*.png"))
    ))

    if roi_overlay_files:
        for overlay_path in roi_overlay_files:
            try:
                img = RLImage(overlay_path)
                max_w, max_h = 175 * mm, 105 * mm
                scale = min(max_w / img.imageWidth, max_h / img.imageHeight, 1.0)
                img.drawWidth = img.imageWidth * scale
                img.drawHeight = img.imageHeight * scale

                caption = os.path.splitext(os.path.basename(overlay_path))[0].replace("_", " ")
                story.append(KeepTogether([
                    Paragraph(caption, h2),
                    img,
                    Spacer(1, 5 * mm)
                ]))
            except Exception as overlay_error:
                story.append(Paragraph(
                    f"ROI overlay could not be embedded: "
                    f"{os.path.basename(overlay_path)} ({overlay_error})",
                    small
                ))
    else:
        story.append(Paragraph(
            "No sphere/ROI overlay PNG files were found at report-generation time.",
            body
        ))

    story.append(PageBreak())

    # ---------- Per-scan QA pages ----------
    story.append(Paragraph("Per-scan QA records", h1))
    story.append(Paragraph(
        "Each record reports the already-calculated features and the Stage III voting result. "
        "The GSR value is merged from the corresponding MI_*.csv file when available.", body
    ))
    story.append(PageBreak())

    voting_path_col = "Pathes" if "Pathes" in voting_df.columns else None

    for modality in ["anat", "diff", "func"]:
        df = feature_tables.get(modality)
        if df is None or df.empty:
            continue
        id_col = _qa_find_id_column(df)
        mi_df = mi_tables.get(modality)
        mi_id_col = _qa_find_id_column(mi_df) if mi_df is not None and not mi_df.empty else None

        # Build a fast lookup for GSR/motion from the auxiliary output.
        mi_lookup = {}
        if mi_df is not None and not mi_df.empty and mi_id_col is not None:
            for _, mi_row in mi_df.iterrows():
                mi_lookup[str(mi_row[mi_id_col])] = mi_row

        for row_number, (_, row) in enumerate(df.iterrows(), start=1):
            scan_id = str(row[id_col])
            votes = pd.DataFrame()
            if voting_path_col is not None:
                votes = voting_df[voting_df[voting_path_col].astype(str) == scan_id]

            if votes.empty:
                vote_count = 0
                stage3_status = "No Stage III outlier vote"
            else:
                vote_count = int(pd.to_numeric(votes.iloc[0].get("Voting outliers (from 5)", 0), errors="coerce") or 0)
                stage3_status = f"Flagged by {vote_count} of 5 methods"

            story.append(Paragraph(f"{modality.upper()} scan {row_number}", scan_title))
            story.append(Paragraph(scan_id.replace("&", "&amp;"), small))
            story.append(Spacer(1, 2 * mm))

            # Select the most useful existing metrics without assuming all modalities have all columns.
            metric_candidates = [
                "SNR Chang", "SNR Normal", "tSNR (Averaged Brain ROI)",
                "Displacement factor (std of Mutual information)", "Ghosting",
                "SpatRx", "SpatRy", "SpatRz"
            ]
            scan_rows = [["Metric", "Value"]]
            for metric in metric_candidates:
                if metric in row.index:
                    scan_rows.append([metric, _qa_numeric(row[metric])])

            mi_row = mi_lookup.get(scan_id)
            if mi_row is not None:
                # Accept the current auxiliary column names without changing upstream naming.
                for col in ["Motion", "Ghosting", "GSR", "Ghosting GSR"]:
                    if col in mi_row.index:
                        label = "GSR (continuous ghosting)" if col in ["Ghosting", "GSR", "Ghosting GSR"] else "Motion (MI std)"
                        # Avoid duplicate Motion/GSR rows if an alias occurs.
                        if not any(r[0] == label for r in scan_rows):
                            scan_rows.append([label, _qa_numeric(mi_row[col])])

            scan_rows.append(["Stage III", stage3_status])
            scan_table = Table(scan_rows, colWidths=[95 * mm, 70 * mm], repeatRows=1)
            scan_table.setStyle(TableStyle([
                ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#D9EAF7")),
                ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
                ("FONTSIZE", (0, 0), (-1, -1), 7.5),
                ("GRID", (0, 0), (-1, -1), 0.35, colors.HexColor("#CCCCCC")),
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
                ("TOPPADDING", (0, 0), (-1, -1), 4),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
            ]))
            if vote_count > 0:
                scan_table.setStyle(TableStyle([
                    ("BACKGROUND", (0, -1), (-1, -1), colors.HexColor("#FCE8E6")),
                    ("FONTNAME", (0, -1), (-1, -1), "Helvetica-Bold")
                ]))
            story.append(scan_table)

            if not votes.empty:
                v = votes.iloc[0]
                vote_cols = ["One_class_SVM", "IsolationForest", "LocalOutlierFactor", " EllipticEnvelope", "statistical_method"]
                vote_rows = [["Stage III indicator", "Flag"]]
                for vc in vote_cols:
                    if vc in v.index:
                        vote_rows.append([vc.strip(), "Yes" if bool(v[vc]) else "No"])
                vt = Table(vote_rows, colWidths=[95 * mm, 70 * mm], repeatRows=1)
                vt.setStyle(TableStyle([
                    ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#EEEEEE")),
                    ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
                    ("FONTSIZE", (0, 0), (-1, -1), 7.5),
                    ("GRID", (0, 0), (-1, -1), 0.35, colors.HexColor("#CCCCCC")),
                    ("TOPPADDING", (0, 0), (-1, -1), 4),
                    ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
                ]))
                story.append(Spacer(1, 3 * mm))
                story.append(vt)

            story.append(PageBreak())

    # ---------- Full Stage III table ----------
    story.append(Paragraph("Stage III outlier table", h1))
    if voting_df.empty:
        story.append(Paragraph("No Stage III flagged scans were present in votings.csv.", body))
    else:
        display_cols = [c for c in [
            "sequence_type", "Pathes", "corresponding_img",
            "One_class_SVM", "IsolationForest", "LocalOutlierFactor",
            " EllipticEnvelope", "statistical_method", "Voting outliers (from 5)"
        ] if c in voting_df.columns]
        tab = [[Paragraph(str(c).strip(), small) for c in display_cols]]
        for _, r in voting_df[display_cols].iterrows():
            tab.append([Paragraph(str(r[c]).replace("&", "&amp;"), small) for c in display_cols])
        available_w = 178 * mm
        col_w = available_w / max(len(display_cols), 1)
        t = Table(tab, colWidths=[col_w] * len(display_cols), repeatRows=1)
        t.setStyle(TableStyle([
            ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#D9EAF7")),
            ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
            ("FONTSIZE", (0, 0), (-1, -1), 5.7),
            ("VALIGN", (0, 0), (-1, -1), "TOP"),
            ("GRID", (0, 0), (-1, -1), 0.25, colors.HexColor("#CCCCCC")),
            ("TOPPADDING", (0, 0), (-1, -1), 2),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 2),
        ]))
        story.append(t)
    story.append(PageBreak())

    # ---------- Provenance ----------
    story.append(Paragraph("Provenance and file integrity", h1))
    story.append(Paragraph(
        "The hashes below identify the tabular inputs used to compose this PDF and support reproducible archiving of the QA protocol.", body
    ))
    provenance_files = feature_files + mi_files + ([voting_path] if os.path.isfile(voting_path) else [])
    prov_rows = [["File", "SHA256"]]
    for f in provenance_files:
        try:
            prov_rows.append([os.path.basename(f), _qa_sha256(f)])
        except Exception:
            prov_rows.append([os.path.basename(f), "hash unavailable"])
    if len(prov_rows) == 1:
        prov_rows.append(["-", "No report input files found"])
    prov_rows_wrapped = [[Paragraph(str(c), small) for c in r] for r in prov_rows]
    prov_table = Table(prov_rows_wrapped, colWidths=[58 * mm, 108 * mm], repeatRows=1)
    prov_table.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#D9EAF7")),
        ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
        ("FONTSIZE", (0, 0), (-1, -1), 6.5),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("GRID", (0, 0), (-1, -1), 0.35, colors.HexColor("#CCCCCC")),
        ("TOPPADDING", (0, 0), (-1, -1), 3),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
    ]))
    story.append(prov_table)

    doc.build(story, onFirstPage=_footer, onLaterPages=_footer)
    return report_path


#%% Adjusting the existing feature table by adding a new sheet to it with the data that need to be discarded

def QCtable(Path, format_type):
    
    ML_algorythms= ML(Path, format_type)
    ML_algorythms=pd.concat(ML_algorythms) 
    ML_algorythms[['One_class_SVM',' EllipticEnvelope','IsolationForest',"LocalOutlierFactor"]]=ML_algorythms[['One_class_SVM',' EllipticEnvelope','IsolationForest',"LocalOutlierFactor"]]==-1 
    Abook = []
    Names =[]
    for file in glob.glob(os.path.join(Path, '*calculated_features*.csv')) :
        
        if "diff" in file:
            dti_path= file
            dti_features= pd.read_csv(dti_path)
            Abook.append(dti_features)
            Names.append("diff")
        elif "func" in file :
            fmri_path= file
            fmri_features= pd.read_csv(fmri_path)
            Abook.append(fmri_features)
            Names.append("func")
        elif "anat" in file :    
             t2w_path= file
             t2w_features= pd.read_csv(t2w_path)
             Abook.append(t2w_features)
             Names.append("anat")    

    
    
  
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
 

         
         

    
 
    
    #prepare outliers
    statiscal=[True if path in Pathes else False for path in ML_algorythms["Pathes"] ]

            
    ML_algorythms["statistical_method"]= statiscal    
    ML_number=list(ML_algorythms[["One_class_SVM" ,'IsolationForest',"LocalOutlierFactor",' EllipticEnvelope',"statistical_method"]].sum(axis=1))
    if format_type == "raw":              
        ML_algorythms= ML_algorythms[["Pathes","sequence_name", "corresponding_img","sequence_type","One_class_SVM" ,'IsolationForest',"LocalOutlierFactor",' EllipticEnvelope',"statistical_method"]]
    elif format_type == "nifti":
        ML_algorythms= ML_algorythms[["Pathes","corresponding_img","sequence_type","One_class_SVM" ,'IsolationForest',"LocalOutlierFactor",' EllipticEnvelope',"statistical_method"]]
    ML_algorythms["Voting outliers (from 5)"]=   ML_number 
    ML_algorythms= ML_algorythms[ML_algorythms["Voting outliers (from 5)"]>=1]
    final_result = os.path.join(Path,"votings.csv")
    ML_algorythms.to_csv(final_result)

    # Stage III is complete. Generate the QA report from the outputs that
    # already exist; report generation does not alter QC calculations.
    try:
        report_path = GenerateQAReport(Path, format_type)
        print("AIDAqc QA report saved to:")
        print(report_path)
    except Exception as report_error:
        # A report-generation problem must not invalidate a completed QC run.
        print("Warning: AIDAqc QA PDF report could not be generated:")
        print(report_error)


    

    

    
 
    
 
#%%  



#%% For Questions please Contact: aref.kalantari-sarcheshmeh@uk-koeln.de




