"""PDF QA report adapted from phantom; reads outputs without changing QC decisions."""
import glob
import hashlib
import os
import subprocess
from datetime import datetime
from html import escape

import pandas as pd

def _qa_sha256(file_path):
    """Return SHA256 for provenance without changing the source file."""
    h = hashlib.sha256()
    with open(file_path, "rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def _qa_git_commit(start_path):
    """Return the current Git commit when the output resides in/under a repository."""
    candidates = [os.path.dirname(os.path.abspath(__file__)), os.path.abspath(start_path)]
    for candidate in candidates:
        try:
            return subprocess.check_output(
                ["git", "-C", candidate, "rev-parse", "--short", "HEAD"],
                stderr=subprocess.DEVNULL,
                universal_newlines=True,
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

    The report reads existing caculated_features_*.csv (or calculated_features_*.csv), MI_*.csv, votings.csv,
    and QCfigures/*.png files. It does not recalculate features or change any
    QC/outlier decisions.
    """
    from reportlab.lib import colors
    from reportlab.lib.enums import TA_CENTER
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
    # Current runs write at the output root; completed runs archive tables below it.
    # Prefer current tables over archived ones, and corrected spelling over legacy spelling.
    feature_by_modality = {}
    for folder in [os.path.join(output_path, "calculated_features"), output_path]:
        for pattern in ["caculated_features_*.csv", "calculated_features_*.csv"]:
            for file_path in sorted(glob.glob(os.path.join(folder, pattern))):
                feature_by_modality[_qa_modality_from_name(file_path)] = file_path
    feature_files = sorted(feature_by_modality.values())
    mi_files = sorted(glob.glob(os.path.join(output_path, "MI_*.csv")))
    voting_path = os.path.join(output_path, "votings.csv")
    if not os.path.isfile(voting_path):
        raise FileNotFoundError("Run Stage III before generating the QA report: votings.csv is missing")
    voting_df = pd.read_csv(voting_path)

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
    run_rows = [[Paragraph(escape(str(c)), small) for c in r] for r in run_rows]
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
    methods_rows_wrapped = [[Paragraph(escape(str(c)), small) for c in r] for r in methods_rows]
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
                    Paragraph(escape(caption), h2),
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
                count = pd.to_numeric(votes.iloc[0].get("Voting outliers (from 5)", 0), errors="coerce")
                vote_count = 0 if pd.isna(count) else int(count)
                stage3_status = f"Flagged by {vote_count} of 5 methods"

            story.append(Paragraph(f"{modality.upper()} scan {row_number}", scan_title))
            story.append(Paragraph(escape(scan_id), small))
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
                for col in ["Motion", "GSR", "Ghosting GSR", "Ghosting"]:
                    if col in mi_row.index:
                        label = "GSR (continuous ghosting)" if col in ["Ghosting", "GSR", "Ghosting GSR"] else "Motion (MI std)"
                        # Avoid duplicate Motion/GSR rows if an alias occurs.
                        if not any(r[0] == label for r in scan_rows):
                            scan_rows.append([label, _qa_numeric(mi_row[col])])

            scan_rows.append(["Stage III", stage3_status])
            scan_rows = [[Paragraph(escape(str(c)), small) for c in r] for r in scan_rows]
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
                        vote_rows.append([vc.strip(), "Yes" if str(v[vc]).strip().lower() in ("true", "1", "1.0") else "No"])
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
        tab = [[Paragraph(escape(str(c).strip()), small) for c in display_cols]]
        for _, r in voting_df[display_cols].iterrows():
            tab.append([Paragraph(escape(str(r[c])), small) for c in display_cols])
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
    prov_rows_wrapped = [[Paragraph(escape(str(c)), small) for c in r] for r in prov_rows]
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
