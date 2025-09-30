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
from sklearn.preprocessing import StandardScaler

import numpy as np
import os
import glob
import warnings
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import MaxNLocator
from scipy import ndimage, signal
import changSNR as ch


# =========================
# Suppress benign scikit-learn warnings
# =========================
warnings.filterwarnings(
    "ignore",
    message="The covariance matrix associated to your dataset is not full rank",
    category=UserWarning,
    module="sklearn.covariance._robust_covariance"
)
warnings.filterwarnings(
    "ignore",
    message="n_neighbors .* is greater than the total number of samples",
    category=UserWarning,
    module="sklearn.neighbors._lof"
)


# =========================
# Compatibility stubs (keep ParsingData.py happy)
# =========================
def tic():
    """Stub for compatibility (no-op timer)."""
    pass

def toc(*args, **kwargs):
    """Stub for compatibility (no-op timer)."""
    pass


# =========================
# Vectorized mutual information
# =========================
def mutualInfo(Im1, Im2, bins=20):
    # Vectorized joint histogram using integer bin indices + bincount
    x = Im1.ravel()
    y = Im2.ravel()

    x_min, x_max = np.min(x), np.max(x)
    y_min, y_max = np.min(y), np.max(y)
    if x_min == x_max:
        x_max = x_min + 1.0
    if y_min == y_max:
        y_max = y_min + 1.0

    x_idx = np.floor((x - x_min) * (bins / (x_max - x_min))).astype(np.int64)
    y_idx = np.floor((y - y_min) * (bins / (y_max - y_min))).astype(np.int64)
    x_idx = np.clip(x_idx, 0, bins - 1)
    y_idx = np.clip(y_idx, 0, bins - 1)

    joint_idx = x_idx * bins + y_idx
    h2d = np.bincount(joint_idx, minlength=bins * bins).reshape(bins, bins)

    pxy = h2d.astype(np.float64)
    s = pxy.sum()
    if s == 0:
        return 0.0
    pxy /= s
    px = pxy.sum(axis=1, keepdims=True)
    py = pxy.sum(axis=0, keepdims=True)
    pxpy = px @ py  # outer product
    nz = pxy > 0
    return np.sum(pxy[nz] * (np.log(pxy[nz]) - np.log(pxpy[nz])))


# =========================
# Ghosting (ou)
# =========================
def GhostCheck(input_file):
    """
    Detect Nyquist ghosting in MR images using mutual information (MI).

    Parameters
    ----------
    input_file : nibabel image
        NIfTI or similar object with get_fdata().

    Returns
    -------
    bool
        True if ghosting detected, False otherwise.
    """
    img_data = input_file.get_fdata()
    if img_data.ndim > 3:
        img_data = img_data.mean(axis=-1)

    H, W, Z = img_data.shape
    im_ref = img_data[:, :, Z // 2]

    # Candidate ghost shifts: W/2, W/4, W/8, ...
    shifts = []
    d = W // 2
    while d >= 1:
        shifts.append(d)
        d //= 2
    shifts = np.unique(np.clip(np.asarray(shifts, dtype=int), 1, W - 1))

    if shifts.size == 0:
        return False

    cols = np.arange(W)
    mi_vals = []
    for s in shifts:
        idx = (cols - s) % W
        shifted = im_ref[:, idx]
        mi_vals.append(mutualInfo(shifted, im_ref))
    mi_vals = np.asarray(mi_vals)

    if mi_vals.size == 0:
        return False

    # --- Detection thresholds ---
    med = np.median(mi_vals)
    mad = np.median(np.abs(mi_vals - med)) + 1e-12  # avoid zero
    height_thr = 0.80 * np.max(mi_vals)   # 80% of max
    prom_thr = 0.0 * mad                  # 0×MAD prominence

    hits = (mi_vals >= height_thr) & ((mi_vals - med) >= prom_thr)

    return bool(np.any(hits))


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

    if img.ndim < 3 or img.mean() == 0:
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

    # Vectorized 8-corner noise regions
    x = int(np.ceil(S[0] * 0.15)); y = int(np.ceil(S[1] * 0.15)); z = int(np.ceil(S[2] * 0.15))
    blocks = (
        img[:x, :y, :z],      img[:x, -y:, :z],      img[-x:, :y, :z],      img[-x:, -y:, :z],
        img[:x, :y, -z:],     img[:x, -y:, -z:],     img[-x:, :y, -z:],     img[-x:, -y:, -z:]
    )
    noise_std = np.std(np.concatenate([b.ravel() for b in blocks]))
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
    if IM.ndim == 3:
        IM = IM[:, :, None, :]

    img = IM.astype("float64")
    fff = 0 if img.shape[-1] < 10 else 10

    sig = img[:, :, :, fff:]
    tsnr_map = 20 * np.log10(sig.mean(axis=-1) / (sig.std(axis=-1) + 1e-12))  # eps for stability

    img_avg = img.mean(axis=-1)
    COM = [int(i) for i in ndimage.center_of_mass(img_avg)]
    r = int(np.floor(0.10 * np.mean(img.shape[:2])))
    Mask = sphere(img.shape[:3], r, COM)
    return np.mean(tsnr_map[Mask])


# =========================
# Motion (rsfMRI)
# =========================
def Ismotion(input_file):
    IM = np.asanyarray(input_file.dataobj)
    if IM.ndim == 3:
        IM = IM[:, :, None, :]

    img = IM.astype("float64")
    fff = 0 if img.shape[-1] < 11 else 10
    sig = img[:, :, :, fff:]

    # Choose slice with highest mean intensity across time
    temp_mean = sig.mean(axis=(0, 1, 3))
    k = int(np.argmax(temp_mean))
    stack = sig[:, :, k, :]  # (H, W, T)
    ref = stack[:, :, 0]

    mi_vals = [mutualInfo(ref, stack[:, :, t]) for t in range(1, stack.shape[-1])]
    final = np.asarray(mi_vals)
    if final.size == 0:
        return final, str([0, 0]), 0.0, 0.0

    max_mov_between = str([final.argmin() + 10, final.argmax() + 10])
    gmv = final.max() - final.min()
    lmv = final.std()
    return final, max_mov_between, gmv, lmv


# =========================
# QC plotting & pies
# =========================
def QCPlot(Path):
    qc_fig_path = os.path.join(Path, "QCfigures")
    if not os.path.isdir(qc_fig_path):
        os.mkdir(qc_fig_path)

    books, names = [], []
    for file in glob.glob(os.path.join(Path, "*caculated_features*.csv")):
        if "diff" in file:
            books.append(pd.read_csv(file)); names.append("diff")
        elif "func" in file:
            books.append(pd.read_csv(file)); names.append("func")
        elif "anat" in file:
            books.append(pd.read_csv(file)); names.append("anat")

    title_font = {"family": "serif", "fontname": "DejaVu Sans"}
    label_font = {"family": "serif", "fontname": "DejaVu Sans"}

    hh = 1
    for df, N in zip(books, names):
        cols = list(df.columns)
        if cols:
            cols.pop(0)

        for C in cols:
            if C in (
                "SNR Chang",
                "tSNR (Averaged Brain ROI)",
                "SNR Normal",
                "Displacement factor (std of Mutual information)",
            ):
                data = pd.to_numeric(df[C], errors="coerce").to_numpy()
                data = data[np.isfinite(data)]
                n = data.size
                if n < 2:
                    continue

                q75, q25 = np.percentile(data, [75, 25])
                iqr = q75 - q25
                rng = data.max() - data.min()

                if rng <= 0 or iqr <= 0:
                    B = 1
                else:
                    h = 2 * iqr / (n ** (1 / 3))
                    B = 1 if (not np.isfinite(h) or h <= 0) else int(np.ceil(rng / h))

                nbins = max(1, B * 7)
                XX = min(22, max(1, B) * 5)

                plt.figure(hh, figsize=(9, 5), dpi=300)
                ax2 = plt.subplot(1, 1, 1, label="hist")
                y, x, bars = plt.hist(data, bins=nbins, histtype="bar", edgecolor="white")

                plt.xlabel(f"{N}: {C} [a.u.]", fontdict=label_font)
                plt.ylabel("Frequency", fontdict=label_font)
                ax2.spines["right"].set_visible(False)
                ax2.spines["top"].set_visible(False)
                plt.locator_params(axis="x", nbins=XX)
                ax2.yaxis.set_major_locator(MaxNLocator(integer=True))

                if iqr > 0:
                    if C == "Displacement factor (std of Mutual information)":
                        ll = q75 + 1.5 * iqr
                        annotate = ("Q3 + 1.5*IQ", lambda bar: bar.get_x() > ll)
                    else:
                        ll = q25 - 1.5 * iqr
                        annotate = ("Q1 - 1.5*IQ", lambda bar: bar.get_x() < ll)

                    if len(y):
                        plt.text(ll, 2 * max(y) / 3, annotate[0], color="grey", fontdict=label_font)
                    for bar in bars:
                        if annotate[1](bar):
                            bar.set_facecolor("red")
                    plt.axvline(ll, color="grey", linestyle="--")

                plt.suptitle(f"{N}: {C}", fontdict=title_font)
                legend = plt.legend(
                    handles=[mpatches.Patch(color="tab:blue", label="Keep"),
                             mpatches.Patch(color="red", label="Discard")],
                    fontsize=8,
                )
                for t in legend.get_texts():
                    t.set_fontfamily("serif"); t.set_fontsize(8)

                ax2.xaxis.set_tick_params(labelsize=8)
                ax2.yaxis.set_tick_params(labelsize=8)

                out = os.path.join(qc_fig_path, f"{C}{N}.png")
                plt.savefig(out, dpi=300); plt.close()

        hh += 1

    # --- Spatial resolution pies ---
    plt.figure(hh, figsize=(9, 5), dpi=300)
    rr = 1
    for df, N in zip(books, names):
        cols = list(df.columns)
        if cols:
            cols.pop(0)
        for C in cols:
            if C in ("SpatRx", "SpatRy", "SpatRz"):
                vals = pd.to_numeric(df[C], errors="coerce").to_numpy()
                vals = vals[np.isfinite(vals)]
                if vals.size == 0:
                    continue
                labels, counts = np.unique(vals, return_counts=True)
                labels2 = [f"{l:.3f} mm" for l in labels]

                ax1 = plt.subplot(len(names), 3, rr)
                ax1.pie(counts, labels=labels2, autopct="%1.0f%%", startangle=180)
                ax1.axis("equal")
                ax1.set_title(f"{N}:{C}", fontdict=title_font)
                plt.suptitle("Resolution homogeneity between data", weight="bold")
                ax1.xaxis.set_tick_params(labelsize=8)
                ax1.yaxis.set_tick_params(labelsize=8)
                rr += 1

    out = os.path.join(qc_fig_path, "Spatial_Resolution.png")
    plt.savefig(out, dpi=300); plt.close()


# =========================
# ML voting
# =========================
def ML(Path, format_type):
    results = []
    for csv in glob.glob(os.path.join(Path, "*_features_*.csv")):
        A = pd.read_csv(csv)
        A = A.dropna(how="all", axis="columns")

        # Build feature matrix
        if format_type == "raw":
            meta_cols = {"path": 1, "seq": 2, "img": 3}
            X = A.iloc[:, 7:].copy()
        else:  # "nifti" or fallback
            meta_cols = {"path": 1, "img": 2}
            X = A.iloc[:, 6:].copy()

        # Coerce to numeric, drop inf/NaN rows and all-NaN cols
        X = X.apply(pd.to_numeric, errors="coerce").replace([np.inf, -np.inf], np.nan)
        X = X.dropna(how="all", axis="columns").dropna(how="any", axis="index")

        # Align metadata
        idx = X.index
        if len(idx) == 0:
            # No usable rows; nothing to add for this file
            continue
        address = A.iloc[idx, meta_cols["path"]].tolist()
        sequence_name = A.iloc[idx, meta_cols["seq"]].tolist() if "seq" in meta_cols else None
        img_name = A.iloc[idx, meta_cols["img"]].tolist() if "img" in meta_cols else None

        # Drop zero-variance columns
        if hasattr(X, "var"):
            X = X.loc[:, X.var(axis=0, ddof=0) > 0]

        n_samples, n_features = X.shape

        # If no usable feature columns remain, default to all inliers and move on
        if n_features == 0:
            df = pd.DataFrame(
                {
                    "One_class_SVM": np.ones(n_samples, dtype=int),
                    " EllipticEnvelope": np.ones(n_samples, dtype=int),
                    "IsolationForest": np.ones(n_samples, dtype=int),
                    "LocalOutlierFactor": np.ones(n_samples, dtype=int),
                }
            )
            if "diff" in csv:
                df["sequence_type"] = "diff"
            elif "func" in csv:
                df["sequence_type"] = "func"
            elif "anat" in csv:
                df["sequence_type"] = "anat"
            df["Pathes"] = address
            if sequence_name is not None:
                df["sequence_name"] = sequence_name
            if img_name is not None:
                df["corresponding_img"] = img_name
            results.append(df)
            continue

        # Scale features for stability
        X = pd.DataFrame(StandardScaler().fit_transform(X), index=X.index)

        # Default predictions = inliers (+1)
        pred_svm = np.ones(n_samples, dtype=int)
        pred_ell = np.ones(n_samples, dtype=int)
        pred_iso = np.ones(n_samples, dtype=int)
        pred_lof = np.ones(n_samples, dtype=int)

        # Fit models with guards
        try:
            if n_samples >= 2:
                pred_svm = OneClassSVM(gamma="auto", kernel="poly", nu=0.05, shrinking=False).fit(X).predict(X)
        except Exception:
            pass
        try:
            if n_samples >= 2:
                pred_ell = EllipticEnvelope(contamination=0.025, random_state=1).fit_predict(X)
        except Exception:
            pass
        try:
            if n_samples >= 2:
                pred_iso = IsolationForest(
                    n_estimators=100, max_samples="auto", contamination=0.05,
                    max_features=1.0, bootstrap=False, n_jobs=-1, random_state=1
                ).fit_predict(X)
        except Exception:
            pass
        try:
            if n_samples >= 2:
                n_neighbors = max(2, min(20, n_samples - 1))
                pred_lof = LocalOutlierFactor(
                    n_neighbors=n_neighbors, algorithm="auto", metric="minkowski",
                    contamination=0.04, novelty=False, n_jobs=-1
                ).fit_predict(X)
        except Exception:
            pass

        df = pd.DataFrame(
            np.vstack([pred_svm, pred_ell, pred_iso, pred_lof]).T,
            columns=["One_class_SVM", " EllipticEnvelope", "IsolationForest", "LocalOutlierFactor"],
        )

        # Metadata
        if "diff" in csv:
            df["sequence_type"] = "diff"
        elif "func" in csv:
            df["sequence_type"] = "func"
        elif "anat" in csv:
            df["sequence_type"] = "anat"
        df["Pathes"] = address
        if sequence_name is not None:
            df["sequence_name"] = sequence_name
        if img_name is not None:
            df["corresponding_img"] = img_name

        results.append(df)

    return results

# =========================
# QC table generation
# =========================
def QCtable(Path, format_type):
    algs = ML(Path, format_type)
    if not algs:
        return
    algs = pd.concat(algs, ignore_index=True)

    # convert -1 flags to boolean outlier flags
    cols_out = ["One_class_SVM", " EllipticEnvelope", "IsolationForest", "LocalOutlierFactor"]
    algs[cols_out] = algs[cols_out] == -1

    books, names = [], []
    for file in glob.glob(os.path.join(Path, "*caculated_features*.csv")):
        if "diff" in file:
            books.append(pd.read_csv(file)); names.append("diff")
        elif "func" in file:
            books.append(pd.read_csv(file)); names.append("func")
        elif "anat" in file:
            books.append(pd.read_csv(file)); names.append("anat")

    pathes = []
    ST, COE, AvV, V, Med, MaX, MiN = [], [], [], [], [], [], []

    for df, N in zip(books, names):
        COL = df.columns
        for C in COL:
            D = pd.to_numeric(df[C], errors="coerce").replace([np.inf, -np.inf], np.nan)

            if C in ("SNR Chang", "tSNR (Averaged Brain ROI)", "SNR Normal"):
                q75, q25 = np.nanpercentile(D, [75, 25])
                iqr = q75 - q25
                Index = (D < (q25 - 1.5 * iqr)) if iqr > 0 else pd.Series(False, index=D.index)

            elif C == "Displacement factor (std of Mutual information)":
                q75, q25 = np.nanpercentile(D, [75, 25])
                iqr = q75 - q25
                Index = (D > (q75 + 1.5 * iqr)) if iqr > 0 else pd.Series(False, index=D.index)

            else:
                continue

            P = df[COL[1]][Index]
            pathes.extend(P)
            ST.extend([N] * len(P))
            COE.extend([C] * len(P))
            AvV.extend([D.mean()] * len(P))
            V.extend(D[Index])
            Med.extend([D.median()] * len(P))
            MiN.extend([D.min()] * len(P))
            MaX.extend([D.max()] * len(P))

    # mark statistical outliers in ML voting
    outlier_set = set(pathes)
    algs["statistical_method"] = algs["Pathes"].isin(outlier_set)
    vote_cols = ["One_class_SVM", "IsolationForest", "LocalOutlierFactor", " EllipticEnvelope", "statistical_method"]
    algs["Voting outliers (from 5)"] = algs[vote_cols].sum(axis=1)
    algs = algs[algs["Voting outliers (from 5)"] >= 1]

    if format_type == "raw":
        keep_cols = ["Pathes", "sequence_name", "corresponding_img", "sequence_type"] + vote_cols
    else:
        keep_cols = ["Pathes", "corresponding_img", "sequence_type"] + vote_cols
    algs = algs[keep_cols]

    final_result = os.path.join(Path, "votings.csv")
    algs.to_csv(final_result, index=False)

#%% For Questions please Contact: aref.kalantari-sarcheshmeh@uk-koeln.de
