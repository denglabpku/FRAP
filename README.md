# FRAP
Fluorescence recovery after photobleaching (FRAP) analysis pipeline
----

## Overview
Fluorescence recovery after photobleaching (FRAP) experiments record the recovery dynamics of fluorescently labelled molecules after photobleached by a strong laser pulse. In this repository, our pipeline is only focusing on the FRAP experiments within live cells. More details about experiment consideration and mathematically modelling of FRAP results have been extensively discussed in many literature, we listed a few at the end of this file for reference.
## Prerequisite
Before you start, make sure you have installed below ImageJ plugins that required for running the ImageJ macro. The package can be directly installed by using Plugins -> Macros -> Install. The purpose of the first two package is for automatic alignment of stacked images. The Bio-Format is to read Nikon's _nd2_ image format. For more details, please go to their websites. For convenience, we provide the installation packages of the first two packages.<br>
1) StackReg-master.zip http://bigwww.epfl.ch/thevenaz/stackreg/
2) turboreg.zip http://bigwww.epfl.ch/thevenaz/turboreg/
3) Bio-Format <br>
Follow the installation instruction if your ImageJ is not pre-installed this package: https://docs.openmicroscopy.org/bio-formats/6.6.0/users/imagej/installing.html

## Steps of FRAP analysis
### 1. Data format conversion (optional)
we usually peform the FRAP experiments on the Nikon A1R confocal microsocpe, the image files are saved as _.nd2_ format. To compatiable with downstream analysis, run the ImageJ macro 'batch_convert_nd2_to_tif.ijm' to convert the file format from _nd2_ to _tif_.
### 2. Align the image stacks
Run the ImageJ macro 'batch_align_stacks_StackReg.ijm' to align image stacks generated from the step 1. This step is to reduce the potential cell drifts.
### 3. Quantify the FRAP intensity from raw images
Run the script 'a1_SegmentQuantifyMovie.m' to analyze the aligned images from step 2, follow the interactive GUI to select the photobleach center, background center, and control center.
### 4. FRAP curve fitting and figure plotting
Run the script 'a2_AnalyzeQuantifiedFRAPdata.m' to combine data from multiple analyzed and normalized FRAP curve from step 3, and to do either one reaction-dominant fitting, or two reaction-dominant fitting, using either linear or logarithmically spaced time axis. To properly choose the fitting model, please see citations for further discussion. This script will also generate figures and fitting results. (Optional) Use the script 'a3_Plot_MultipleFRAPcurves.m' to compare FRAP curves from multiple samples/conditions, and visualize them on one figure.

## Citation
1) McNally, James G. “Quantitative FRAP in Analysis of Molecular Binding Dynamics In Vivo.” In Methods in Cell Biology, 85:329–51. Elsevier, 2008. https://doi.org/10.1016/S0091-679X(08)85014-5.
2) Sprague, Brian L., Robert L. Pego, Diana A. Stavreva, and James G. McNally. “Analysis of Binding Reactions by Fluorescence Recovery after Photobleaching.” Biophysical Journal 86, no. 6 (June 2004): 3473–95. https://doi.org/10.1529/biophysj.103.026765.
3) Hansen, Anders S., Iryna Pustova, Claudia Cattoglio, Robert Tjian, and Xavier Darzacq. “CTCF and Cohesin Regulate Chromatin Loop Stability with Distinct Dynamics.” ELife 6 (May 3, 2017): e25776. https://doi.org/10.7554/eLife.25776.

