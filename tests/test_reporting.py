"""Run with: python -m unittest discover -s tests -v (requirements-dev.txt)."""
import contextlib
import hashlib
import io
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

os.environ.setdefault("MPLBACKEND", "Agg")
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import nibabel as nib
import numpy as np
import pandas as pd
import pymupdf
from PIL import Image

import FeatureCheck as fc
import QC
from qa_report import GenerateQAReport

MOTION = "Displacement factor (std of Mutual information)"
VOTES = ["One_class_SVM", " EllipticEnvelope", "IsolationForest", "LocalOutlierFactor", "statistical_method"]


def image(ghost=0):
    data = np.ones((40, 40, 8), dtype=float)
    data[4:36, 0:16, :] = np.linspace(100, 200, 512).reshape(32, 16, 1)
    data[4:36, 20:36, :] += ghost
    return nib.Nifti1Image(data, np.eye(4))


class GhostingTests(unittest.TestCase):
    def test_gsr_increases_with_half_fov_replica(self):
        self.assertAlmostEqual(QC.GhostGSRValue(image()), 0.0)
        self.assertGreater(QC.GhostGSRValue(image(20)), QC.GhostGSRValue(image(5)))

    def test_4d_uses_temporal_mean_and_does_not_mutate_input(self):
        data = image(10).get_fdata()
        original = data.copy()
        series = nib.Nifti1Image(np.stack([data, 2 * data], axis=-1), np.eye(4))
        mean = nib.Nifti1Image(1.5 * data, np.eye(4))
        self.assertAlmostEqual(QC.GhostGSRValue(series), QC.GhostGSRValue(mean))
        np.testing.assert_array_equal(data, original)

    def test_empty_foreground_is_unavailable(self):
        blank = nib.Nifti1Image(np.zeros((40, 40, 8)), np.eye(4))
        self.assertTrue(np.isnan(QC.GhostGSRValue(blank)))


class FeatureOutputTests(unittest.TestCase):
    def test_nifti_and_bruker_outputs_keep_existing_features(self):
        for raw in [False, True]:
            for modality in ["anat", "diff", "func"]:
                with self.subTest(raw=raw, modality=modality), tempfile.TemporaryDirectory() as tmp:
                    root = Path(tmp)
                    scan = root / "study" / "1"
                    scan.mkdir(parents=True)
                    im = image(10)
                    if raw:
                        (scan / "pdata" / "1").mkdir(parents=True)
                        (scan / "pdata" / "1" / "visu_pars").touch()
                        (scan / "acqp").touch()
                        address = str(scan)
                    else:
                        scan = scan / "scan.nii.gz"
                        nib.save(im, scan)
                        address = str(scan)
                    pd.DataFrame([address]).to_csv(root / ("raw_data_addreses_" + modality + ".csv"), index=False)
                    with contextlib.ExitStack() as stack:
                        stack.enter_context(contextlib.redirect_stdout(io.StringIO()))
                        stack.enter_context(patch.object(fc, "plot_img"))
                        stack.enter_context(patch.object(fc, "GhostCheck", return_value=True))
                        stack.enter_context(patch.object(fc, "GhostGSRValue", return_value=0.123))
                        stack.enter_context(patch.object(fc, "snrCalclualtor_chang", return_value=25.0))
                        stack.enter_context(patch.object(fc, "snrCalclualtor_normal", return_value=30.0))
                        stack.enter_context(patch.object(fc, "TsnrCalclualtor", return_value=12.0))
                        stack.enter_context(patch.object(fc, "Ismotion", return_value=(np.array([0.1, 0.2]), "[0,1]", 0.1, 0.05)))
                        if raw:
                            factory = stack.enter_context(patch.object(fc.pr, "Bruker2Nifti"))
                            factory.return_value.nim = im
                            stack.enter_context(patch.object(fc.par, "read_param_file", return_value=({}, {"ACQ_method": "RARE", "ACQ_protocol_name": "T2W"})))
                            fc.CheckingRawFeatures(str(root))
                        else:
                            fc.CheckingNiftiFeatures(str(root))
                    features = pd.read_csv(root / ("caculated_features_" + modality + ".csv"))
                    aux = pd.read_csv(root / ("MI_" + modality + ".csv"))
                    self.assertTrue(features.loc[0, "Ghosting"])
                    self.assertEqual(aux.loc[0, "Ghosting"], 0.123)
                    self.assertNotIn("GSR", features.columns)
                    self.assertEqual(features.loc[0, "FileAddress"], aux.loc[0, "FileAddress"])
                    if modality == "anat":
                        self.assertTrue(pd.isna(aux.loc[0, "Motion"]))
                        self.assertNotIn(MOTION, features)
                    else:
                        self.assertEqual(features.loc[0, MOTION], aux.loc[0, "Motion"])


class ReportTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)

    def fixtures(self, empty_votes=False):
        rows = []
        for modality in ["anat", "diff", "func"]:
            scan_id = "/data/subject_<A&B>/" + modality + "/" + "long_name_" * 8 + ".nii.gz"
            frame = pd.DataFrame({"FileAddress": [scan_id], "Ghosting": [True], "SpatRx": [0.1], "SpatRy": [0.1], "SpatRz": [0.5], "SNR Chang": [25.0]})
            frame.to_csv(self.root / ("caculated_features_" + modality + ".csv"))
            pd.DataFrame({"FileAddress": [scan_id], "Motion": [np.nan if modality == "anat" else 0.15], "Ghosting": [0.123]}).to_csv(self.root / ("MI_" + modality + ".csv"))
            row = dict(Pathes=scan_id, sequence_type=modality, corresponding_img="scan_<A&B>.png")
            row.update({name: name == "One_class_SVM" for name in VOTES})
            row["Voting outliers (from 5)"] = 1
            rows.append(row)
        voting = pd.DataFrame(rows)
        if empty_votes:
            voting = voting.iloc[:0]
        voting.to_csv(self.root / "votings.csv")
        (self.root / "QCfigures").mkdir()
        Image.new("RGB", (320, 150), "white").save(self.root / "QCfigures" / "QC_<A&B>.png")

    def report_text(self):
        before = {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in self.root.rglob("*.csv")}
        path = GenerateQAReport(str(self.root), "nifti")
        self.assertTrue(Path(path).is_file())
        self.assertEqual(before, {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in before})
        with pymupdf.open(path) as doc:
            self.assertGreaterEqual(len(doc), 6)
            return "\n".join(page.get_text() for page in doc)

    def test_report_contains_all_modalities_metrics_votes_and_hashes(self):
        self.fixtures()
        text = self.report_text()
        for value in ["ANAT scan 1", "DIFF scan 1", "FUNC scan 1", "0.1230", "Flagged by 1 of 5 methods", "SHA256"]:
            self.assertIn(value, text)
        self.assertIn("<A&B>", text)
        self.assertNotIn("ROI placement overlays", text)

    def test_report_handles_empty_vote_table(self):
        self.fixtures(empty_votes=True)
        text = self.report_text()
        self.assertIn("No Stage III outlier vote", text)
        self.assertIn("No Stage III flagged scans", text)

    def test_report_reads_archived_corrected_tables(self):
        self.fixtures()
        archived = self.root / "calculated_features"
        archived.mkdir()
        for p in self.root.glob("caculated_features_*.csv"):
            p.rename(archived / p.name.replace("caculated", "calculated"))
        text = self.report_text()
        self.assertIn("ANAT scan 1", text)
        self.assertIn("FUNC scan 1", text)

    def test_current_table_takes_precedence_over_archived_table(self):
        self.fixtures()
        archived = self.root / "calculated_features"
        archived.mkdir()
        pd.DataFrame({"FileAddress": ["stale_scan"]}).to_csv(archived / "calculated_features_anat.csv")
        self.assertNotIn("stale_scan", self.report_text())

    def test_report_requires_stage_three_results(self):
        with self.assertRaisesRegex(FileNotFoundError, "Stage III"):
            GenerateQAReport(str(self.root), "nifti")

    def test_pdf_failure_preserves_voting_csv(self):
        self.fixtures()
        feature = pd.read_csv(self.root / "caculated_features_anat.csv")
        votes = pd.DataFrame({"Pathes": feature["FileAddress"], "sequence_type": ["anat"], "corresponding_img": ["scan.png"]})
        for col in VOTES[:-1]:
            votes[col] = -1
        with patch.object(QC, "ML", return_value=[votes]), patch.object(QC, "GenerateQAReport", side_effect=RuntimeError("PDF unavailable")), contextlib.redirect_stdout(io.StringIO()):
            QC.QCtable(str(self.root), "nifti")
        result = pd.read_csv(self.root / "votings.csv")
        self.assertEqual(result.loc[0, "Voting outliers (from 5)"], 4)


if __name__ == "__main__":
    unittest.main()
