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
import time
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import MaxNLocator
from scipy import ndimage, signal
from matplotlib import font_manager as fm  # kept in case of external usage
from matplotlib.font_manager import FontProperties  # kept for compatibility
import changSNR as ch


# =========================
# Simple tic/toc utilities
# =========================
def TicTocGenerator():
    ti = 0
    tf = time.time()
    while True:
        ti = tf
        tf = time.time()
        yield tf - ti


TicToc = TicTocGenerator()


def toc(tempBool=True):
    dt = next(TicToc)
    if tempBool:
        print("Elapsed time: %f seconds.\n" % dt)


def tic():
    toc(False)


# =========================
# Ghosting
# =========================
def GhostCheck(input_file):
    img = input_file
    img_data = img.get_fdata()
    img_shape = img_data.shape

    Mmos = []
    n = 1
    while (img_shape[1] % (2 ** n)) == 0:
        Mmos.append(img_shape[1] / 2 ** n)
        n += 1
    Mmos = np.asarray(Mmos)

    if img_data.ndim > 3:
        img_data = img_data.mean(axis=-1)

    Im_ref = img_data[:, :, int(img_shape[2] / 2)]
    MI_vec = []
    for ii in range(int(img_shape[1])):
        Im_rol = np.roll(Im_ref, ii, axis=1)
        MI_vec.append(mutualInfo(Im_rol, Im_ref))

    peaks_strong, _ = signal.find_peaks(MI_vec, height=0.25 * max(MI_vec))
    peaks_weak, _ = signal.find_peaks(MI_vec)

    StrongGhost = np.sum(np.isin(peaks_strong, Mmos))
    WeekGhost = np.sum(np.isin(peaks_weak, Mmos))

    return (WeekGhost > 2) or (StrongGhost > 0)


# =========================
# Resolution
# =========================
def ResCalculator(input_file):
    HDR = input_file.header
    return HDR["pixdim"][1:4]


# =========================
# SNR (Chang)
# =========================
def snrCalclualtor_chang(input_file):
    imgData = input_file
    IM = np.asanyarray(imgData.dataobj)
    imgData = np.squeeze(IM.astype("float64"))

    mm = imgData.mean()
    if mm == 0:
        return np.nan

    if imgData.ndim < 3:
        return np.nan

    snr_chang_slice_vec = []
    ns = imgData.shape[2]
    if imgData.ndim > 3:
        n_dir = imgData.shape[-1]
        if n_dir < 10:
            fff = 0
            print("\nWarning: 4th dim < 10 (unstable values possible)\n")
        else:
            fff = 5

    if ns > 4:
        ns_lower = int(np.floor(ns / 2) - 2)
        ns_upper = int(np.floor(ns / 2) + 2)
    else:
        ns_lower = 0
        ns_upper = ns

    if imgData.ndim == 3:
        for slc in range(ns_lower, ns_upper):
            Slice = imgData[:, :, slc]
            try:
                _, estStdChang, _ = ch.calcSNR(Slice, 0, 1)
            except ValueError:
                estStdChang = np.nan
            snr_chang_slice = 20 * np.log10(np.mean(Slice) / estStdChang)
            snr_chang_slice_vec.append(snr_chang_slice)
    else:
        n_dir = imgData.shape[-1]
        for slc in range(ns_lower, ns_upper):
            for bb in range(fff, n_dir - 1):
                Slice = imgData[:, :, slc, bb]
                try:
                    _, estStdChang, _ = ch.calcSNR(Slice, 0, 1)
                except ValueError:
                    estStdChang = np.nan
                snr_chang_slice = 20 * np.log10(np.mean(Slice) / estStdChang)
                snr_chang_slice_vec.append(snr_chang_slice)

    snr_chang_slice_vec = np.asarray(snr_chang_slice_vec)
    m = np.mean(snr_chang_slice_vec[~np.isinf(snr_chang_slice_vec) & ~np.isnan(snr_chang_slice_vec)])
    return m


# =========================
# SNR (Normal)
# =========================
def snrCalclualtor_normal(input_file):
    IM = np.asanyarray(input_file.dataobj)
    imgData = np.squeeze(IM.astype("float64"))

    if imgData.ndim < 3:
        imgData = np.tile(imgData[:, :, np.newaxis], (1, 1, 10))

    Data = imgData
    S = np.squeeze(Data).shape

    if len(S) == 3:
        imgData = np.squeeze(Data)
    elif len(S) == 4:
        imgData = np.squeeze(Data[:, :, :, 0])

    S = np.squeeze(imgData).shape

    # CoM + spherical mask
    COM = [int(i) for i in ndimage.center_of_mass(imgData)]
    r = int(np.floor(0.10 * np.mean(S)))
    if r > S[2]:
        r = S[2]
    Mask = sphere(S, r, COM)
    Signal = imgData[Mask].mean()

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
    noise_blocks = [np.squeeze(imgData[idx]) for idx in corners]
    Noise_std = np.std(np.array(noise_blocks))

    SNR = 20 * np.log10(Signal / Noise_std)
    if np.isinf(SNR):
        SNR = np.nan
        print("Impossible: SNR = inf (likely zeros around ROI). Replacing with NaN.")
    return SNR


def show_slices(slices):
    fig, axes = plt.subplots(1, len(slices))
    for i, Slice in enumerate(slices):
        axes[i].imshow(Slice.T, cmap="gray", origin="lower")


def sphere(shape, radius, position):
    assert len(position) == len(shape)
    semisizes = (radius,) * len(shape)
    grid = [slice(-x0, dim - x0) for x0, dim in zip(position, shape)]
    position = np.ogrid[grid]
    arr = np.zeros(shape, dtype=float)
    for x_i, semisize in zip(position, semisizes):
        arr += (x_i / semisize) ** 2
    return arr <= 1.0


# =========================
# tSNR
# =========================
def TsnrCalclualtor(input_file):
    IM = np.asanyarray(input_file.dataobj)
    S = IM.shape
    if len(S) == 3:
        IM = IM.reshape((S[0], S[1], 1, S[2]))

    imgData = IM.astype("float64")

    if IM.shape[-1] < 10:
        fff = 0
    else:
        fff = 10

    signal_averge_over_time = imgData[:, :, :, fff:].mean(axis=-1)
    signal_std_over_time = imgData[:, :, :, fff:].std(axis=-1)
    tSNR_map = 20 * np.log10(signal_averge_over_time / signal_std_over_time)

    S = IM.shape
    imgData_average = imgData.mean(axis=-1)

    COM = [int(i) for i in ndimage.center_of_mass(imgData_average)]
    r = int(np.floor(0.10 * np.mean([S[0:2]])))
    Mask = sphere(S[0:3], r, COM)
    tSNR = np.mean(tSNR_map[Mask])

    return tSNR


# =========================
# Mutual Information
# =========================
def mutualInfo(Im1, Im2):
    t1_slice = Im1
    t2_slice = Im2

    hist_2d, x_edges, y_edges = np.histogram2d(t1_slice.ravel(), t2_slice.ravel(), bins=20)

    hist_2d_log = np.zeros_like(hist_2d)
    non_zeros = hist_2d != 0
    hist_2d_log[non_zeros] = np.log(hist_2d[non_zeros])

    pxy = hist_2d / float(np.sum(hist_2d))
    px = np.sum(pxy, axis=1)
    py = np.sum(pxy, axis=0)
    px_py = px[:, None] * py[None, :]
    nzs = pxy > 0
    MI = np.sum(pxy[nzs] * np.log(pxy[nzs] / px_py[nzs]))
    return MI


# =========================
# Motion detection (rsfMRI)
# =========================
def Ismotion(input_file):
    imgData = input_file
    IM = np.asanyarray(imgData.dataobj)
    S = IM.shape
    if len(S) == 3:
        IM = IM.reshape((S[0], S[1], 1, S[2]))

    if IM.shape[-1] < 11:
        fff = 0
    else:
        fff = 10

    imgData = IM[:, :, :, fff:].astype("float64")
    S = imgData.shape
    temp_mean = imgData.mean(axis=(0, 1, 3))
    temp_max = temp_mean.argmax()
    temp_Data = imgData[:, :, temp_max, :]
    Im_fix = temp_Data[:, :, 0]
    Im_rot = temp_Data

    MI_all = []
    for z in range(1, S[-1]):
        MI_all.append(mutualInfo(Im_fix, Im_rot[:, :, z]))

    Final = np.asarray(MI_all)
    Max_mov_between = str([Final.argmin() + 10, Final.argmax() + 10])
    GMV = getrange(Final)
    LMV = np.std(Final)

    return Final, Max_mov_between, GMV, LMV


def getrange(numbers):
    return max(numbers) - min(numbers)


# =========================
# Plotting QC histograms & resolution pies
# =========================
def QCPlot(Path):
    saving_path = Path
    QC_fig_path = os.path.join(Path, "QCfigures")
    if not os.path.isdir(QC_fig_path):
        os.mkdir(QC_fig_path)

    Abook = []
    Names = []
    for file in glob.glob(os.path.join(Path, "*caculated_features*.csv")):
        if "diff" in file:
            Abook.append(pd.read_csv(file))
            Names.append("diff")
        elif "func" in file:
            Abook.append(pd.read_csv(file))
            Names.append("func")
        elif "anat" in file:
            Abook.append(pd.read_csv(file))
            Names.append("anat")

    hh = 1
    rr = 1
    title_font = {"family": "serif", "fontname": "DejaVu Sans"}
    label_font = {"family": "serif", "fontname": "DejaVu Sans"}

    for nn, N in enumerate(Names):
        COL = list(Abook[nn].columns)
        if COL:
            COL.pop(0)
        D = Abook[nn]

        for cc, C in enumerate(COL):
            if C in (
                "SNR Chang",
                "tSNR (Averaged Brain ROI)",
                "SNR Normal",
                "Displacement factor (std of Mutual information)",
            ):
                # --- histogram ---
                plt.figure(hh, figsize=(9, 5), dpi=300)
                ax2 = plt.subplot(1, 1, 1, label="histogram")

                # Coerce numeric + finite
                Data = pd.to_numeric(D[C], errors="coerce").to_numpy()
                Data = Data[np.isfinite(Data)]
                if Data.size == 0:
                    plt.close()
                    continue

                q75, q25 = np.percentile(Data, [75, 25])
                iqr = q75 - q25
                data_range = Data.max() - Data.min()

                # Freedman–Diaconis with guards; ensure >=1 bin
                if data_range <= 0 or iqr <= 0:
                    B = 1
                else:
                    h = 2 * iqr / (Data.size ** (1 / 3))
                    if not np.isfinite(h) or h <= 0:
                        B = 1
                    else:
                        B = int(np.ceil(data_range / h))

                XX = min(22, max(1, B) * 5)
                y, x, bars = plt.hist(
                    Data, bins=max(1, B * 7), histtype="bar", edgecolor="white"
                )

                plt.xlabel(f"{N}: {C} [a.u.]", fontdict=label_font)
                plt.ylabel("Frequency", fontdict=label_font)
                ax2.spines["right"].set_visible(False)
                ax2.spines["top"].set_visible(False)
                plt.locator_params(axis="x", nbins=XX)
                ax2.yaxis.set_major_locator(MaxNLocator(integer=True))

                # Outlier threshold only if iqr > 0
                if iqr > 0:
                    if C == "Displacement factor (std of Mutual information)":
                        ll = q75 + 1.5 * iqr
                        plt.text(
                            1.07 * ll,
                            2 * max(y) / 3,
                            "Q3 + 1.5*IQ",
                            color="grey",
                            fontdict=label_font,
                        )
                        for bar in bars:
                            if bar.get_x() > ll:
                                bar.set_facecolor("red")
                    else:
                        ll = q25 - 1.5 * iqr
                        plt.text(
                            1.001 * ll,
                            2 * max(y) / 3,
                            "Q1 - 1.5*IQ",
                            color="grey",
                            fontdict=label_font,
                        )
                        for bar in bars:
                            if bar.get_x() < ll:
                                bar.set_facecolor("red")
                    plt.axvline(ll, color="grey", linestyle="--")

                plt.suptitle(f"{N}: {C}", fontdict=title_font)

                red_patch = mpatches.Patch(color="red", label="Discard")
                blue_patch = mpatches.Patch(color="tab:blue", label="Keep")
                legend = plt.legend(handles=[blue_patch, red_patch], fontsize=8)
                for text in legend.get_texts():
                    text.set_fontfamily("serif")
                    text.set_fontsize(8)

                ax2.xaxis.set_tick_params(labelsize=8)
                ax2.yaxis.set_tick_params(labelsize=8)

                base_filename = os.path.join(QC_fig_path, C + N)
                plt.savefig(base_filename + ".png", dpi=300)
                plt.close()

        hh += 1

    # --- Spatial resolution pies ---
    plt.figure(hh, figsize=(9, 5), dpi=300)
    for nn, N in enumerate(Names):
        COL = list(Abook[nn].columns)
        if COL:
            COL.pop(0)
        D = Abook[nn]
        for cc, C in enumerate(COL):
            if C in ("SpatRx", "SpatRy", "SpatRz"):
                Data = pd.to_numeric(D[C], errors="coerce").to_numpy()
                Data = Data[np.isfinite(Data)]
                if Data.size == 0:
                    continue

                labels = np.unique(Data)
                sizes = [(Data == l).sum() for l in labels]
                labels_rounded = np.round(labels, 3)
                labels2 = [f"{l} mm" for l in labels_rounded]

                ax1 = plt.subplot(len(Names), 3, rr)
                ax1.pie(sizes, labels=labels2, autopct="%1.0f%%", startangle=180)
                ax1.axis("equal")
                ax1.set_title(f"{N}:{C}", fontdict=title_font)
                plt.suptitle("Resolution homogeneity between data", weight="bold")
                ax1.xaxis.set_tick_params(labelsize=8)
                ax1.yaxis.set_tick_params(labelsize=8)
                rr += 1

    base_filename = os.path.join(QC_fig_path, "Spatial_Resolution")
    plt.savefig(base_filename + ".png", dpi=300)
    plt.close()


# =========================
# ML voting
# =========================
def ML(Path, format_type):
    result = []
    for N, csv in enumerate(glob.glob(os.path.join(Path, "*_features_*.csv"))):
        csv_path = os.path.join(Path, csv)
        Abook = pd.read_csv(csv_path)

        if np.any(Abook.isnull().all()[:]):
            print("The following csv file contains NaN values for one or more of its features:")
            print(csv_path)
            print("Voting can not be conducted.")
            print("Analyzing next sequence...")
            continue

        Abook = Abook.dropna(how="all", axis="columns")
        Abook = Abook.dropna(how="any")

        address = [i for i in Abook.iloc[:, 1]]

        if format_type == "raw":
            sequence_name = [i for i in Abook.iloc[:, 2]]
            img_name = [i for i in Abook.iloc[:, 3]]
            X = Abook.iloc[:, 7:]
        elif format_type == "nifti":
            img_name = [i for i in Abook.iloc[:, 2]]
            X = Abook.iloc[:, 6:]
        else:
            # default fallback
            img_name = [i for i in Abook.iloc[:, 2]]
            X = Abook.iloc[:, 6:]

        # Models
        nu = 0.05
        clf = OneClassSVM(gamma="auto", kernel="poly", nu=nu, shrinking=False).fit(X)
        svm_pre = clf.predict(X)

        elpenv = EllipticEnvelope(contamination=0.025, random_state=1)
        ell_pred = elpenv.fit_predict(X)

        iforest = IsolationForest(
            n_estimators=100,
            max_samples="auto",
            contamination=0.05,
            max_features=1.0,
            bootstrap=False,
            n_jobs=-1,
            random_state=1,
        )
        iso_pred = iforest.fit_predict(X)

        lof = LocalOutlierFactor(
            n_neighbors=20,
            algorithm="auto",
            metric="minkowski",
            contamination=0.04,
            novelty=False,
            n_jobs=-1,
        )
        local_pred = lof.fit_predict(X)

        algorythms = [svm_pre, ell_pred, iso_pred, local_pred]
        arr = np.vstack(algorythms).T
        df = pd.DataFrame(
            arr,
            columns=["One_class_SVM", " EllipticEnvelope", "IsolationForest", "LocalOutlierFactor"],
        )

        if "diff" in csv:
            df["sequence_type"] = ["diff"] * len(df)
        elif "func" in csv:
            df["sequence_type"] = ["func"] * len(df)
        elif "anat" in csv:
            df["sequence_type"] = ["anat"] * len(df)

        df["Pathes"] = address
        if format_type == "raw":
            df["sequence_name"] = sequence_name
        df["corresponding_img"] = img_name

        result.append(df)

    return result


def QCtable(Path, format_type):
    ML_algorythms = ML(Path, format_type)
    if not ML_algorythms:
        return
    ML_algorythms = pd.concat(ML_algorythms)

    # convert -1 flags to boolean outlier flags
    cols_out = ["One_class_SVM", " EllipticEnvelope", "IsolationForest", "LocalOutlierFactor"]
    ML_algorythms[cols_out] = ML_algorythms[cols_out] == -1

    Abook = []
    Names = []
    for file in glob.glob(os.path.join(Path, "*caculated_features*.csv")):
        if "diff" in file:
            Abook.append(pd.read_csv(file))
            Names.append("diff")
        elif "func" in file:
            Abook.append(pd.read_csv(file))
            Names.append("func")
        elif "anat" in file:
            Abook.append(pd.read_csv(file))
            Names.append("anat")

    ST, COE, AvV, V, Pathes, Med, MaX, MiN = [], [], [], [], [], [], [], []

    for nn, N in enumerate(Names):
        d = Abook[nn]
        COL = d.columns

        for cc, C in enumerate(COL):
            D = d[C]

            if C in ("SNR Chang", "tSNR (Averaged Brain ROI)", "SNR Normal"):
                D = pd.to_numeric(D, errors="coerce")
                D = D.replace([np.inf, -np.inf], np.nan)
                q75, q25 = np.nanpercentile(D, [75, 25])
                iqr = q75 - q25
                if iqr > 0:
                    ll = q25 - 1.5 * iqr
                    Index = D < ll
                else:
                    Index = pd.Series([False] * len(D), index=D.index)

                P = d[COL[1]][Index]
                M, Me, Mi, Ma = D.mean(), D.median(), D.min(), D.max()

                Pathes.extend(P)
                ST.extend([N] * len(P))
                COE.extend([C] * len(P))
                AvV.extend([M] * len(P))
                V.extend(D[Index])
                Med.extend([Me] * len(P))
                MiN.extend([Mi] * len(P))
                MaX.extend([Ma] * len(P))

            if C == "Displacement factor (std of Mutual information)":
                D = pd.to_numeric(D, errors="coerce")
                D = D.replace([np.inf, -np.inf], np.nan)
                q75, q25 = np.nanpercentile(D, [75, 25])
                iqr = q75 - q25
                if iqr > 0:
                    ul = q75 + 1.5 * iqr
                    Index = D > ul
                else:
                    Index = pd.Series([False] * len(D), index=D.index)

                P = d[COL[1]][Index]
                M, Me, Mi, Ma = D.mean(), D.median(), D.min(), D.max()

                Pathes.extend(P)
                ST.extend([N] * len(P))
                COE.extend([C] * len(P))
                AvV.extend([M] * len(P))
                V.extend(D[Index])
                Med.extend([Me] * len(P))
                MiN.extend([Mi] * len(P))
                MaX.extend([Ma] * len(P))

            if N == "ErrorData":
                D = d[C]
                Pathes.extend(D)
                S = "Faulty Data"
                ST.extend([S] * len(D))
                COE.extend(["-"] * len(D))
                AvV.extend(["-"] * len(D))
                V.extend(["-"] * len(D))
                Med.extend(["-"] * len(D))
                MiN.extend(["-"] * len(D))
                MaX.extend(["-"] * len(D))

    # prepare outliers for voting
    statiscal = [path in set(Pathes) for path in ML_algorythms["Pathes"]]
    ML_algorythms["statistical_method"] = statiscal

    ML_number = list(ML_algorythms[["One_class_SVM", "IsolationForest", "LocalOutlierFactor", " EllipticEnvelope", "statistical_method"]].sum(axis=1))

    if format_type == "raw":
        ML_algorythms = ML_algorythms[
            [
                "Pathes",
                "sequence_name",
                "corresponding_img",
                "sequence_type",
                "One_class_SVM",
                "IsolationForest",
                "LocalOutlierFactor",
                " EllipticEnvelope",
                "statistical_method",
            ]
        ]
    elif format_type == "nifti":
        ML_algorythms = ML_algorythms[
            [
                "Pathes",
                "corresponding_img",
                "sequence_type",
                "One_class_SVM",
                "IsolationForest",
                "LocalOutlierFactor",
                " EllipticEnvelope",
                "statistical_method",
            ]
        ]

    ML_algorythms["Voting outliers (from 5)"] = ML_number
    ML_algorythms = ML_algorythms[ML_algorythms["Voting outliers (from 5)"] >= 1]

    final_result = os.path.join(Path, "votings.csv")
    ML_algorythms.to_csv(final_result, index=False)

#%% For Questions please Contact: aref.kalantari-sarcheshmeh@uk-koeln.de




