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

See the full manual [here](https://github.com/Aswendt-Lab/AIDAqc/blob/main/docs/AIDAqc_v2_1.pdf).

For installation in a [apptainer](https://apptainer.org/) container for GNU/Linux:

```{bash}

# Download the repository
git clone https://github.com/Aswendt-Lab/AIDAqc.git
cd AIDAqc

# Create a new apptainer container
apptainer build aidaqc.sif apptainer.def

# Get into a bash shell in the container
apptainer run aidaqc.sif

```

<h3>The story behind this tool</h3> 

It can be challenging to acquire MR images of consistent quality or to decide between good vs. bad quality data in large databases. Manual screening without quantitative criteria is strictly user-dependent and for large databases is neither practical nor in the spirit of good scientific practice. In contrast to clinical MRI, in animal MRI, there is no consensus on the standardization of quality control measures or categorization of good vs. bad quality images. As we were forced to screen hundreds of scans for a recent project, we decided to automate this process as part of our Atlas-based Processing Pipeline (AIDA).

<h3>Validation and Datasets</h3> 

This tool has been validated and used in the following publication: [Publication Link](https://gin.g-node.org/Aswendt_Lab/2023_Kalantari_AIDAqc)

A total of 23 datasets from various institutes were used for validation and testing. These datasets can be found via: [Datasets Link](https://gin.g-node.org/Aswendt_Lab/2023_Kalantari_AIDAqc)

<h3>Download test dataset</h3>
https://gin.g-node.org/Aswendt_Lab/testdata_aidaqc

[<h3><b>CONTACT</h3></b>](https://neurologie.uk-koeln.de/forschung/ag-neuroimaging-neuroengineering/)
Aref Kalantari (aref.kalantari-sarcheshmehATuk-koeln.de) and Markus Aswendt (markus.aswendtATuk-koeln.de)

<h3><b>LICENSE</h3></b>

[GNU General Public License v3.0](https://github.com/aswendtlab/AIDAqc/blob/main/LICENSE)
