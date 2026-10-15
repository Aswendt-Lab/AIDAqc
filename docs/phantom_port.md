# Phantom improvements transferred to open-dev

Comparison baseline: main `be03a1d`, open-dev `0e87305`, phantom `581cd34`.

| Area | main / original open-dev | phantom | open-dev transfer |
| --- | --- | --- | --- |
| SNR and tSNR ROI | Animal-MRI sphere based on 10% of mean volume dimensions | 20% of mean in-plane dimensions; thin-stack ellipsoid fallback | Preserve animal-MRI calculations |
| Motion | Standard deviation of mutual information for diffusion and functional scans | Also exported in MI CSVs | Export the existing values |
| Ghosting | Boolean peak-based flag | Additional continuous intensity-based GSR | Add GSR in auxiliary CSVs only |
| QA PDF | Absent | Report after Stage III | Port into qa_report.py, invoked after voting |
| ROI overlay images | Absent | Sphere, noise ROI, and ellipsoid overlays | Keep phantom-specific overlays on phantom |
| Feature filenames | Legacy caculated_features spelling | Corrected calculated_features spelling | Preserve existing outputs; report accepts both |
| Documentation and contacts | open-dev lacked recent main updates | Phantom-specific manual | Bring forward main Markdown manual and contact updates |

The report preserves the existing voting decisions and records the tabular inputs with SHA256 hashes. Its generator is separate from QC calculations. It uses QCfigures PNGs when available and summarizes per-scan results; phantom ROI-overlay sections are omitted because this workflow does not generate those images.

Auxiliary MI CSVs retain the phantom schema: FileAddress, Motion, Ghosting (plus the saved CSV index). Ghosting in these files is continuous GSR, not the boolean Ghosting column in the main feature tables. Anatomical Motion is NaN. No new column is added to the machine-learning feature matrix.

GSR is an image-coordinate estimate using a single central slice and a half-FOV shift along image axis 1. It does not infer phase encoding from acquisition metadata and should be interpreted after verifying image orientation. The transfer does not establish new animal-MRI validation or change any exclusion threshold.

The legacy environment and Dockerfile are retained. No phantom Conda environment or phantom manual is substituted for the animal-MRI setup.

The percentile foreground mask can include ghost signal when the actual object occupies a small fraction of the field of view, which can suppress the reported GSR. Inspect the images alongside the score; a zero GSR does not establish absence of ghosting. The transferred calculation is unchanged from phantom.

Regression checks: install requirements-dev.txt and run `python -m unittest discover -s tests -v`. The suite covers both raw and NIfTI auxiliary output wiring, controlled GSR sensitivity and temporal averaging, report content and input-file integrity, legacy/corrected archived tables, and preservation of voting output after report failure. The output-wiring tests isolate the established metric functions; they do not validate animal-MRI quality thresholds.
