# JAKSTATModels
https://doi.org/10.5281/zenodo.14545524

# Modeling of JAK/STAT signaling activity in response to IL4 stimulation

## Introduction
This repository contains code and data files used in the preparation of the manuscript "Capturing the dynamics of STAT6 macrophage polarization using bioluminescence temporal spectrometry" by Zheng et al.

The following authors contributed to this repository:

Belinda S Akpa (University of Tennessee and National Institute for Modeling Biological Systems)

Zachary Fox (Oak Ridge National Laboratory)

### Code contents
JAKSTATModel_InferredFromBTSWayne.m

### Data file contents

data_CtrlAndIL4Replicates.csv

Figure7BParameterValues.mat

## Usage
The code file, JAKSTATModel_InferredFromBTSWayne.m, produces the plot shown in Figure 7B of the relevant manuscript. The code uses the best parameter set identified in the ABC-SMC search to generate response curves for IL4 stimulation.  These curves are then plotted against the experimental data.  

### JAKSTATModel_InferredFromBTSWayne.m
Contains a function JAKSTATModel_InferredFromBTSWayne(parameterValues) with the following argument:
- parameterValues - 1x26 array containing the estimated parameter values for the model, expressed as log10(theta). For Figure 7B, the corresponding parameter values are found in the file Figure7BParameterValues.mat.

