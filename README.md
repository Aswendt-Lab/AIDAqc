<img align="left" src="https://github.com/maswendt/AIDAmri/blob/master/AIDA_Logo.png" width="120">
<h1>AIDA<i>qc</i></h1>

*An automated and simple tool for fast qualtiy analysis of animal MRI*
<br/>
<br/>
<h3>Features</h3> 
<img align="center" src="https://github.com/aswendtlab/AIDAqc/blob/main/AIDAqc_workflow.pdf>
- **Input:** Bruker raw data or NIFTY (T2-weighted MRI, diffusion weighted MRI or DTI, and rs-fMRI)
- **Calculations:** SNR, tSNR, movement variability, data quality categorization (finds bad quality outlier)
- **Output Format:** Excel sheets & pdf

[**See the poster for all details**](https://github.com/aswendtlab/AIDAqc/blob/main/AIDAqc_Poster_Summary.pdf) 

<h3>Installation</h3> 
Download the repository => Install Python 3.6 (Anaconda) => Import AIDAqc conda environment aidaqc.yaml

Two main functions: *ParsingAllrawData* and *CheckingFeatures_final*

See the full manual [here](https://github.com/aswendtlab/AIDAqc/blob/main/AIDAqc_help_v1_1.pdf).

<h3>The story behind this tool</h3> 

It can be challenging to acquire MR images of consistent quality or to decide for the good vs. bad qualtiy data in a large databases. Manual screening without quantitative criteria is strictly user-dependent and for large databases is neither practical nor in the spirit of good scientific practice. In contrast to clinical MRI, in animal MRI, there is no consensus on standardization of quality control measures or categorization of good vs. bad quality images. As we were forced for a recent project to sreen hundreds of scans, we decided to automate this processa as part of our Atlas-based Processing Pipeline (AIDA).

<h3>Download test dataset</h3>
https://gin.g-node.org/arefks/AIDAqc_test_data

[<h3><b>CONTACT</h3></b>](https://neurologie.uk-koeln.de/forschung/ag-neuroimaging-neuroengineering/)
Aref Kalantari (arefks@gmail.com) and Markus Aswendt (markus.aswendt@uk-koeln.de)

<h3><b>LICENSE</h3></b>

[GNU General Public License v3.0](https://github.com/aswendtlab/AIDAqc/blob/main/LICENSE)
