# Wilson's guide to 10x Multiome GEX + ATAC analysis  
Date: 12/2/2021  
Author: Wilson  

## Description  
I use R package Seurat to run most of the analysis here.
This tutorial assumes that you have generated the output from cellranger arc.
If not, I will create additional tutorial to morph other possible input for the code.
Right now, the tutorial assumes the input files per sample to be:  
- matrix h5 (filtered_feature_bc_matrix.h5)  
- fragment file (fragments.tsv.gz)  


## Project 1: Differential GEX + ATAC analysis between two conditions  
Loading required packages:  

```
library(Signac) ## needed for single-cell chromatin analysis
library(Seurat) ## base package for single-cell analysis 
library(ggplot2) ## for plotting
library(JASPAR2020) ## needed for ChromVAR analysis
library(TFBSTools) ## needed for ChromVAR analysis
library(EnsDb.Rnorvegicus.v79) ## needed for one function (dont have to be exact, unless you are running chromvar)
library(BSgenome.Rnorvegicus.UCSC.rn6) ## needed for one function (dont have to be exacr, unless you are running chromvar)
```  


