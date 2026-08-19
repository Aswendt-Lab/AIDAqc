<img align="right" src="https://github.com/Aswendt-Lab/AIDAqc/blob/main/docs/AIDA_Logo_wide.001.png" width="500">
<h1>AIDA<i>qc</i></h1>

*An automated and simple tool for fast quality analysis of animal MRI*
<br/>
<br/>
<h3>Features</h3> 

- **Input:** Bruker raw data or NIFTY (T2-weighted MRI, diffusion-weighted MRI, or DTI, and rs-fMRI)
- **Calculations:** SNR, tSNR, movement variability, data quality categorization (finds bad quality outliers)
- **Output Format:** CSV sheets, PDFs, & images

<img align="left" src="https://github.com/Aswendt-Lab/AIDAqc/blob/main/docs/AIDAqc_workflow.png">

<br/>
<br/>

[**See the poster for all details**](https://github.com/Aswendt-Lab/AIDAqc/blob/main/docs/AIDAqc_Poster_Summary.pdf) 

<h3>Installation</h3> 
Download the repository => Install Python 3.6 (Anaconda) => Import AIDAqc conda environment aidaqc.yaml

Main function: *ParsingData*

See the full manual [here](https://github.com/Aswendt-Lab/AIDAqc/blob/main/docs/AIDAqc_v2_2.pdf).

<h3>Docker/Apptainer Usage</h3>

AIDAqc provides a multi-architecture Docker image supporting AMD64 
(Intel/AMD Linux) and ARM64 (e.g., Apple Silicon).

<details>
<summary><b>Build the Docker image</b></summary>

```bash
docker buildx build \
  --platform linux/amd64,linux/arm64 \
  -t aswendtlab/aidaqc:2.2 \
  -t aswendtlab/aidaqc:latest \
  --push .
```

The appropriate Conda environment is selected automatically:

- `linux/amd64` → `aidaqc-intel.yaml`
- `linux/arm64` → `aidaqc-arm64.yaml`

</details>

<details>
<summary><b>Check available options</b></summary>

```bash
docker run --rm aswendtlab/aidaqc:2.2 -h
```

</details>

<details>
<summary><b>Run AIDAqc with NIfTI data</b></summary>

```bash
docker run --rm \
  -v /your/project/data:/data \
  -v /your/project/qc:/qc \
  aswendtlab/aidaqc:2.2 \
  -i /data \
  -o /qc \
  -f nifti
```

</details>

<details>
<summary><b>Run AIDAqc with Bruker raw data</b></summary>

```bash
docker run --rm \
  -v /your/project/data:/data \
  -v /your/project/qc:/qc \
  aswendtlab/aidaqc:2.2 \
  -i /data \
  -o /qc \
  -f raw
```

</details>

All generated AIDAqc outputs, including calculated feature tables, 
outlier-detection results, QC figures, ROI-placement figures, and 
`AIDAqc_QA_Report.pdf`, are written to the mounted `/qc` directory.

<details>
<summary><b>Apptainer</b></summary>

```bash
git clone https://github.com/Aswendt-Lab/AIDAqc.git
cd AIDAqc

apptainer build aidaqc.sif apptainer.def
apptainer shell aidaqc.sif
```

</details>

<h3>Tutorial</h3>

The [YouTube tutorial](https://youtu.be/SP4sWW313DQ?si=4WaTI544FzAkBVbY) guides you through the workflow (note: this is for v1.0).

<h3>The story behind this tool</h3> 

It can be challenging to acquire MR images of consistent quality or to decide between good vs. bad quality data in large databases. Manual screening without quantitative criteria is strictly user-dependent and for large databases is neither practical nor in the spirit of good scientific practice. In contrast to clinical MRI, in animal MRI, there is no consensus on the standardization of quality control measures or categorization of good vs. bad quality images. As we were forced to screen hundreds of scans for a recent project, we decided to automate this process as part of our Atlas-based Processing Pipeline (AIDA).

<h3>Validation and Datasets</h3> 

This tool has been validated and used in the following publication: [Publication Link](https://gin.g-node.org/Aswendt_Lab/2023_Kalantari_AIDAqc)

A total of 23 datasets from various institutes were used for validation and testing. These datasets can be found via: [Datasets Link](https://gin.g-node.org/Aswendt_Lab/2023_Kalantari_AIDAqc)

<h3>Download test dataset</h3>

[Dataset Link](https://gin.g-node.org/Aswendt_Lab/testdata_aida)

If you encounter problems, report directly in [![Gitter](https://badges.gitter.im/AIDA_tools/community.svg)](https://gitter.im/AIDA_tools/community?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge)
or 
join our Open Office Hour - each Thursday 3:00 pm (UTC+2) [![Zoom](https://img.shields.io/badge/Zoom-2D8CFF?style=for-the-badge&logo=zoom&logoColor=white)](https://uni-frankfurt.zoom-x.de/j/63112745009?pwd=JBTjMVbuaTw9cZvFnppTwCPjGdQEyx.1)


For all other inquiries: Markus Aswendt (aswendtATmed.uni-frankfurt.de)

<h3><b>LICENSE</h3></b>

[GNU General Public License v3.0](https://github.com/aswendtlab/AIDAqc/blob/main/LICENSE)
