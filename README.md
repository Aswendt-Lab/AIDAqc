<h1>AIDA<i>qc</i></h1>
<img align="right" src="https://github.com/Aswendt-Lab/AIDAqc/blob/main/docs/AIDA_Logo_wide.001.png" width="500">

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

See the full manual [here](docs/AIDAqc_v2_1.md).

<h3>Docker/Apptainer Usage</h3>

```{bash}
#Build

docker build aidaqc:2.1 .

# Running the main ParsingData.py:

docker run -v /your/project/data:/data -v /your/project/qc aidaqc:2.1 -i /data -o /qc -f raw

```

For installation in a [apptainer](https://apptainer.org/) container for GNU/Linux:
```{bash}
# Download the repository
git clone https://github.com/Aswendt-Lab/AIDAqc.git
cd AIDAqc

# Create a new apptainer container
apptainer build aidaqc.sif apptainer.def

# Get into a bash shell in the container
apptainer shell aidaqc.sif

```
<h3>Branches</h3>

AIDAqc provides separate workflows for animal MRI and MRI phantom data:

- [**`main`**](https://github.com/Aswendt-Lab/AIDAqc/tree/main) – the stable animal MRI workflow.
- [**`phantom`**](https://github.com/Aswendt-Lab/AIDAqc/tree/phantom) – specifically for **MRI phantom data only**. The former `open-dev` implementation has moved here, including phantom-specific ROI analysis, additional motion and ghosting outputs, and the PDF QA report. See the [phantom README](https://github.com/Aswendt-Lab/AIDAqc/blob/phantom/README.md) for setup and usage, and the [updated manual](https://github.com/Aswendt-Lab/AIDAqc/blob/phantom/docs/AIDAqc_open_dev.pdf) for the workflow details.
- [**`open-dev`**](https://github.com/Aswendt-Lab/AIDAqc/tree/open-dev) – general development, aligned with `main` after moving the phantom workflow. PDF export will be ported to this baseline in a later change; the phantom-specific analysis stays in `phantom`.

Researchers and developers are welcome to fork the repository and submit pull requests against the branch appropriate to their workflow.

<h3>Tutorial</h3>

To guide you through running the pipeline, please watch the [YouTube tutorial](https://youtu.be/SP4sWW313DQ?si=4WaTI544FzAkBVbY).

<h3>The story behind this tool</h3> 

It can be challenging to acquire MR images of consistent quality or to decide between good vs. bad quality data in large databases. Manual screening without quantitative criteria is strictly user-dependent and for large databases is neither practical nor in the spirit of good scientific practice. In contrast to clinical MRI, in animal MRI, there is no consensus on the standardization of quality control measures or categorization of good vs. bad quality images. As we were forced to screen hundreds of scans for a recent project, we decided to automate this process as part of our Atlas-based Processing Pipeline (AIDA).

<h3>Validation and Datasets</h3> 

This tool has been validated and used in the following publication: [Publication Link](https://gin.g-node.org/Aswendt_Lab/2023_Kalantari_AIDAqc)

A total of 23 datasets from various institutes were used for validation and testing. These datasets can be found via: [Datasets Link](https://gin.g-node.org/Aswendt_Lab/2023_Kalantari_AIDAqc)

<h3>Download test dataset</h3>

[Dataset Link](https://gin.g-node.org/Aswendt_Lab/testdata_aida)

<h3>Contact</h3>

Please report bugs and technical issues directly in [GitHub Issues](https://github.com/Aswendt-Lab/AIDAqc/issues).

For all other requests, email [aswendtlab@gmail.com](mailto:aswendtlab@gmail.com).

<h3><b>LICENSE</h3></b>

[GNU General Public License v3.0](https://github.com/aswendtlab/AIDAqc/blob/main/LICENSE)
