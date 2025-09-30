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
import glob
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import MaxNLocator
from scipy import ndimage, signal
import changSNR as ch


# =========================
# Compatibility stubs
# =========================
def tic():
    """Stub for compatibility (no-op timer)."""
    return None

def toc(*args, **kwargs):
    """Stub for compatibility (no-op timer)."""
    return None


# =========================
# Ghosting
# =========================
def GhostCheck(input_file):
    img_data = input_file.get_fdata()
    img_shape = img_data.shape

    # candidate ghost positions (powers of two along phase-encode)
    mmos = []
    n = 1
    while (img_shape[1] % (2 ** n)) == 0:
        mmos.append(img_shape[1] / 2 ** n)
        n += 1
    mmos = np.asarray(mmos)

    if img_data.ndim > 3:
        img_data = img_data.mean(axis=-1)

    im_ref = img_data[:, :, img_shape[2] // 2]
    mi_vec = [mutualInfo(np.roll(im_ref, ii, axis=1), im_ref) for ii in range(img_shape[1])]

    peaks_strong, _ = signal.find_peaks(mi_vec, height=0.25 * np.max(mi_vec))
    peaks_weak, _ = signal.find_peaks(mi_vec)

    strong = np.sum(np.isin(peaks_strong, mmos))
    weak = np.sum(np.isin(peaks_weak, mmos))
    return (weak > 2) or (strong > 0)


# =========================
# Resolution
# =========================
def ResCalculator(input_file):
    return input_file.header["pixdim"][1:4]


# =========================
# SNR (Chang)
# =========================
def snrCalclualtor_chang(input_file):
    IM = np.asanyarray(input_file.dataobj)
    img = np.squeeze(IM.astype("float64"))

    if img.mean() == 0 or img.ndim < 3:
        return np.nan

    ns = img.shape[2]
    if img.ndim > 3:
        n_dir = img.shape[-1]
        fff = 0 if n_dir < 10 else 5

    if ns > 4:
        sl_lo = int(np.floor(ns / 2) - 2)
        sl_hi = int(np.floor(ns / 2) + 2)
    else:
        sl_lo, sl_hi = 0, ns

    vals = []
    if img.ndim == 3:
        for sl in range(sl_lo, sl_hi):
            slc = img[:, :, sl]
            try:
                _, est_std, _ = ch.calcSNR(slc, 0, 1)
            except ValueError:
                est_std = np.nan
            vals.append(20 * np.log10(np.mean(slc) / est_std))
    else:
        n_dir = img.shape[-1]
        for sl in range(sl_lo, sl_hi):
            for bb in range(fff, n_dir - 1):
                slc = img[:, :, sl, bb]
                try:
                    _, est_std, _ = ch.calcSNR(slc, 0, 1)
                except ValueError:
                    est_std = np.nan
                vals.append(20 * np.log10(np.mean(slc) / est_std))

    vals = np.asarray(vals)
    mask = np.isfinite(vals)
    return np.mean(vals[mask]) if mask.any() else np.nan


# =========================
# SNR (Normal)
# =========================
def snrCalclualtor_normal(input_file):
    IM = np.asanyarray(input_file.dataobj)
    img = np.squeeze(IM.astype("float64"))

    if img.ndim < 3:
        img = np.tile(img[:, :, np.newaxis], (1, 1, 10))

    Data = img
    S = np.squeeze(Data).shape

    if len(S) == 3:
        img = np.squeeze(Data)
    elif len(S) == 4:
        img = np.squeeze(Data[:, :, :, 0])

    S = np.squeeze(img).shape

    # Spherical brain mask around center-of-mass
    COM = [int(i) for i in ndimage.center_of_mass(img)]
    r = int(np.floor(0.10 * np.mean(S)))
    r = min(r, S[2])
    Mask = sphere(S, r, COM)
    signal_val = img[Mask].mean()

    x = int(np.ceil(S[0] * 0.15))
    y = int(np.ceil(S[1] * 0.15))
    z = int(np.ceil(S[2] * 0.15))

    corners = [
        (slice(None, x), slice(None, y), slice(None, z)),
        (slice(None, x), slice(-y, None), slice(None, z)),
        (slice(-x, None), slice(None, y), slice(None, z)),
        (slice(-x, None), slice(-y, None), slice(None, z)),
        (slice(None, x), slice(None, y), slice(-z, None)),
        (slice(None, x), slice(-y, None), slice(-z, None)),
        (slice(-x, None), slice(None, y), slice(-z, None)),
        (slice(-x, None), slice(-y, None), slice(-z, None)),
    ]
    noise_blocks = [np.squeeze(img[idx]) for idx in corners]
    noise_std = np.std(np.array(noise_blocks))

    snr = 20 * np.log10(signal_val / noise_std)
    return snr if np.isfinite(snr) else np.nan


def sphere(shape, radius, position):
    """Generate an n-D spherical mask."""
    assert len(position) == len(shape)
    semisizes = (radius,) * len(shape)
    grid = [slice(-x0, dim - x0) for x0, dim in zip(position, shape)]
    axes = np.ogrid[grid]
    arr = np.zeros(shape, dtype=float)
    for ax, semi in zip(axes, semisizes):
        arr += (ax / semi) ** 2
    return arr <= 1.0


# =========================
# tSNR
# =========================
def TsnrCalclualtor(input_file):
    IM = np.asanyarray(input_file.dataobj)
    S = IM.shape
    if len(S) == 3:
        IM = IM.reshape((S[0], S[1], 1, S[2]))

    img = IM.astype("float64")
    fff = 0 if IM.shape[-1] < 10 else 10

    sig_mean = img[:, :, :, fff:].mean(axis=-1)
    sig_std = img[:, :, :, fff:].std(axis=-1)
    tsnr_map = 20 * np.log10(sig_mean / sig_std)

    img_avg = img.mean(axis=-1)
    COM = [int(i) for i in ndimage.center_of_mass(img_avg)]
    r = int(np.floor(0.10 * np.mean([IM.shape[0:2]])))
    Mask = sphere(IM.shape[0:3], r, COM)
    return np.mean(tsnr_map[Mask])


# =========================
# Mutual Information / motion
# =========================
def mutualInfo(Im1, Im2):
    h2d, _, _ = np.histogram2d(Im1.ravel(), Im2.ravel(), bins=20)
    pxy = h2d / float(np.sum(h2d))
    px = pxy.sum(axis=1)
    py = pxy.sum(axis=0)
    px_py = px[:, None] * py[None, :]
    nz = pxy > 0
    return np.sum(pxy[nz] * np.log(pxy[nz] / px_py[nz]))


def Ismotion(input_file):
    IM = np.asanyarray(input_file.dataobj)
    S = IM.shape
    if len(S) == 3:
        IM = IM.reshape((S[0], S[1], 1, S[2]))

    fff = 0 if IM.shape[-1] < 11 else 10
    img = IM[:, :, :, fff:].astype("float64")
    S = img.shape

    temp_mean = img.mean(axis=(0, 1, 3))
    temp_max = temp_mean.argmax()
    stack = img[:, :, temp_max, :]
    im_fix = stack[:, :, 0]

    mi_vals = [mutualInfo(im_fix, stack[:, :, z]) for z in range(1, S[-1])]
    final = np.asarray(mi_vals)
    max_mov_between = str([final.argmin() + 10, final.argmax() + 10])
    gmv = final.max() - final.min()
    lmv = final.std()
    return final, max_mov_between, gmv, lmv

#%% For Questions please Contact: aref.kalantari-sarcheshmeh@uk-koeln.de




